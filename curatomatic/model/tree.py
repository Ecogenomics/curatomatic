import sys
from collections import defaultdict, deque
from pathlib import Path

import dendropy

from curatomatic.model.metadata import Metadata
from curatomatic.model.node_label import NodeLabel
from curatomatic.model.rank import RANKS
from curatomatic.model.red import RED
from curatomatic.model.taxonomy import Taxonomy
from curatomatic.util.log import logger
from curatomatic.util.tree import parse_node_label, get_tax_strings_for_leaf_nodes, transfer_taxon


class Tree:

    def __init__(self, tree: dendropy.Tree):
        self.tree = tree

        # Calculate required attributes
        self.depth_to_nodes: dict[int, set[dendropy.Node]] = self.get_depth_to_nodes(tree)
        self.node_to_descs: dict[dendropy.Node, set[dendropy.Node]] = self.get_node_desc_nodes(self.depth_to_nodes)
        self.node_labels: dict[dendropy.Node, NodeLabel] = self.get_node_labels(tree)
        self.gid_to_node: dict[str, dendropy.Node] = self.get_gid_to_leaf_node(self.node_labels)
        self.node_to_taxon_labels: dict[dendropy.Node, Taxonomy] = self.get_node_to_taxon_labels(self.node_labels)
        self.taxon_labels_to_node: dict[str, dendropy.Node] = self.get_taxon_labels_to_node(self.node_to_taxon_labels)

    def reload(self):
        self.node_labels: dict[dendropy.Node, NodeLabel] = self.get_node_labels(self.tree)
        self.node_to_taxon_labels: dict[dendropy.Node, Taxonomy] = self.get_node_to_taxon_labels(self.node_labels)
        self.taxon_labels_to_node: dict[str, dendropy.Node] = self.get_taxon_labels_to_node(self.node_to_taxon_labels)

    @staticmethod
    def get_node_labels(tree: dendropy.Tree):
        out = dict()
        for node in tree.postorder_node_iter():
            out[node] = NodeLabel.from_node(node)
        return out

    @staticmethod
    def get_node_to_taxon_labels(node_labels: dict[dendropy.Node, NodeLabel]):
        """Just return those with a taxon label."""
        return {k: v.taxonomy for k, v in node_labels.items() if v.taxonomy}

    @staticmethod
    def get_taxon_labels_to_node(node_to_taxon_labels: dict[dendropy.Node, Taxonomy]):
        out = dict()
        for node, taxonomy in node_to_taxon_labels.items():
            for cur_taxon in taxonomy.iter_ranks():
                if cur_taxon in out:
                    logger.error(f'Taxon label {cur_taxon} is not unique on the tree.')
                    sys.exit(1)
                out[cur_taxon] = node
        return out

    @staticmethod
    def get_depth_to_nodes(tree: dendropy.Tree):
        out = defaultdict(set)
        queue = deque([(tree.seed_node, 0)])
        while len(queue) > 0:
            cur_node, cur_depth = queue.popleft()
            out[cur_depth].add(cur_node)
            for child_node in cur_node.child_nodes():
                queue.append((child_node, cur_depth + 1))
        return out

    @staticmethod
    def get_gid_to_leaf_node(node_labels: dict[dendropy.Node, NodeLabel]):
        out = dict()
        for node in node_labels:
            if node.is_leaf():
                gid = node_labels[node].gid
                if gid.startswith('D-'):
                    continue
                assert gid is not None
                out[gid] = node
        return out

    @staticmethod
    def get_node_desc_nodes(depth_to_nodes: dict[int, set[dendropy.Node]]):
        out = defaultdict(set)
        max_depth = max(depth_to_nodes.keys())
        for i in reversed(range(max_depth + 1)):
            for cur_node in depth_to_nodes[i]:
                if cur_node.is_leaf():
                    out[cur_node] = set()
                else:
                    for child_node in cur_node.child_node_iter():
                        out[cur_node].update(out[child_node])
                        out[cur_node].add(child_node)
        return out

    @classmethod
    def from_file(cls, path: Path):
        tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)

        # Convert any placeholds (e.g. g__) that exist on the basal node
        # to an unresolved name
        # TODO: This is a hack, fix it.
        for leaf_node in tree.leaf_node_iter():
            if leaf_node.taxon.label.startswith('D-'):
                continue

            parent_node = leaf_node.parent_node
            red, bs, taxon = parse_node_label(parent_node)
            assert red == 1.0
            assert bs is None

            for cur_taxon in taxon.split(';'):
                cur_taxon = cur_taxon.strip()
                if cur_taxon != 's__':
                    continue
                if len(cur_taxon) == 3:
                    gid = leaf_node.taxon.label.split('|RED')[0]
                    new_label = 's__' + '{unresolved} ' + gid
                    new_label = f'{new_label}|RED={red:.3f}'
                    parent_node.label = new_label

        return cls(tree=tree)

    def to_file(self, path: Path):
        self.tree.write(path=str(path), schema="newick", suppress_rooting=True)

    def n_nodes_and_zombies(self) -> tuple[int, int]:
        n_nodes, n_zombies = 0, 0
        for leaf_node in self.tree.leaf_node_iter():
            if leaf_node.taxon.label.startswith('D-'):
                n_zombies += 1
            else:
                n_nodes += 1
        return n_nodes, n_zombies

    def find_taxa_for_refinement(self, poly, min_bootstrap: float, red_dict: RED):
        """
        Find all taxa that can be moved to a higher node in the tree without
        creating any nested ranks.
        """
        d_taxon_to_candidates = dict()

        # Iterate over each named node
        for node, taxonomy in self.node_to_taxon_labels.items():

            # Process each of the taxa contained in the string
            for cur_taxon in taxonomy.iter_ranks():

                # Not going to bother refining nested groups
                if cur_taxon in poly:
                    continue

                cur_rank = cur_taxon[0]

                # Do not refine species
                if cur_rank == 's':
                    continue

                cur_rank_idx = RANKS.index(cur_rank)

                # Go up the tree to see if this will conflict
                candidates = list()

                # Keep going up the tree until we find a conflict or hit the root
                cur_node = node.parent_node
                while cur_node is not None:

                    # Find all named taxa that are below the current node
                    desc_taxa = defaultdict(set)
                    for desc_node in self.node_to_descs[cur_node]:
                        cur_desc_taxonomy = self.node_labels[desc_node].taxonomy
                        if cur_desc_taxonomy:
                            for cur_desc_taxon in cur_desc_taxonomy.iter_ranks():
                                desc_taxa[cur_desc_taxon[0]].add(cur_desc_taxon)

                    # Would any of these cause a conflict?
                    conflicts = desc_taxa[cur_rank] - {cur_taxon}
                    if len(conflicts) == 0:
                        # Check bootstrap support
                        cur_node_label = self.node_labels[cur_node]
                        bs_supported = cur_node_label.bootstrap >= min_bootstrap

                        # Only consider this node if we are less than the half way
                        # point between the expected RED for the parent rank
                        # and the current rank
                        parent_rank_red = red_dict.data[RANKS[cur_rank_idx - 1]]
                        cur_rank_red = red_dict.data[cur_rank]
                        halfway_red = abs(cur_rank_red - parent_rank_red) / 2 + parent_rank_red

                        red_supported = cur_node_label.red >= halfway_red

                        # Add to candidate list if this is supported
                        if red_supported and bs_supported:
                            candidates.append(cur_node)
                        cur_node = cur_node.parent_node
                    else:
                        # This would create a nested rank, stop going up
                        cur_node = None

                # If we found any candidates for this taxon, then save them
                if len(candidates) > 0:
                    d_taxon_to_candidates[cur_taxon] = candidates
        return d_taxon_to_candidates

    def get_nested_ranks(self):
        """
        This method goes through the tree to find cases where nodes are nested.
        """
        poly = set()
        for leaf_node in self.tree.leaf_node_iter():
            d_seen = defaultdict(list)
            cur_node = leaf_node

            while cur_node is not None:
                cur_node_tax = self.node_labels[cur_node].taxonomy
                if cur_node_tax:
                    for taxon in cur_node_tax.iter_ranks():
                        rank = taxon[0]
                        if len(d_seen[rank]) > 0:
                            # print(d_seen, taxon)
                            poly.add(taxon)
                            poly.update(d_seen[rank])
                        d_seen[rank].append(taxon)
                cur_node = cur_node.parent_node
        return poly

    def get_leaf_nodes_missing_rank(self):
        out = dict()

        for leaf_node in self.tree.leaf_node_iter():

            # Skip zombie nodes
            if leaf_node.taxon.label.startswith('D-'):
                continue

            # Go up the tree and read the tax string from the tree
            taxonomy = Taxonomy()
            cur_node = leaf_node.parent_node
            while cur_node is not None:

                # Get the tax string from this node (if it exists)
                cur_taxonomy = NodeLabel.from_node(cur_node).taxonomy
                if cur_taxonomy:
                    for cur_taxon in cur_taxonomy.iter_ranks():
                        taxonomy.add_taxon(cur_taxon)

                if taxonomy.all_assigned():
                    # Check if we can exit early
                    cur_node = None
                else:
                    # Keep going up
                    cur_node = cur_node.parent_node

            # Check if we obtained a full tax string
            if not taxonomy.all_assigned():
                missing_ranks = taxonomy.missing_ranks()
                out[leaf_node] = missing_ranks
        return out

    def remove_red_labels(self):
        for node in self.tree.postorder_node_iter():
            if node.is_leaf():
                node.taxon.label = node.taxon.label.split('|')[0]
            else:
                node.label = node.label.split('|')[0]
        return

    def find_low_bootstrap_support_nodes_with_taxon_placement(self, min_bs: float, red_dict: RED, metadata: Metadata):
        nodes_to_refine = list()
        for node, taxonomy in self.node_to_taxon_labels.items():
            node_label = self.node_labels[node]

            # Skip those that are not set (i.e. leaf nodes).
            if node_label.bootstrap is None:
                continue

            if node_label.bootstrap < min_bs:
                nodes_to_refine.append(node)

        # Go over each node and find all candidate nodes that it can be moved to
        # not considering any cases like nesting
        node_to_candidates = defaultdict(set)
        for node in nodes_to_refine:
            node_taxonomy = self.node_to_taxon_labels[node]
            node_label = self.node_labels[node]

            for taxon in node_taxonomy.iter_ranks():
                # Get the prefixes of the taxon / parent / child
                cur_rank = taxon[0]

                # Get the median RED dict values for those
                cur_rank_red_target = red_dict.data[cur_rank]

                # Find the delta
                red_delta_to_current = abs(cur_rank_red_target - node_label.red)

                # Check the children
                for desc_node in self.node_to_descs[node]:
                    desc_node_label = self.node_labels[desc_node]

                    # No point checking if the bootstrap is low
                    if desc_node_label.bootstrap is None or desc_node_label.bootstrap < min_bs:
                        continue

                    # Make sure that the original genome that provides this taxon is down this path
                    desc_node_red_delta = abs(desc_node_label.red - cur_rank_red_target)

                    # Check if the child node is better
                    if desc_node_red_delta < red_delta_to_current:
                        node_to_candidates[(node, taxon)].add(desc_node)

                # Check the parents
                cur_node = node.parent_node
                while cur_node is not None:
                    cur_node_label = self.node_labels[cur_node]

                    # No point checking if the bootstrap is low
                    if cur_node_label.bootstrap is None or cur_node_label.bootstrap < min_bs:
                        break

                    cur_node_red_delta = abs(cur_node_label.red - cur_rank_red_target)

                    # Check if the parent node is better
                    if cur_node_red_delta < red_delta_to_current:
                        node_to_candidates[(node, taxon)].add(cur_node)
                        cur_node = cur_node.parent_node
                    else:
                        break


        return

    def refine_existing_labels(self, rank: str, min_bs: float, red_dict: RED, metadata: Metadata):
        """
        Refine an existing rank in the tree, designed to be called form phylum to genus.
        """
        out = defaultdict(list)

        d_gid_to_taxonomy = metadata.get_taxonomy_from_gids(set(self.gid_to_node.keys()))
        d_taxon_to_gids = metadata.get_taxon_to_gids(set(self.gid_to_node.keys()))
        singleton_taxa = {k for k, v in d_taxon_to_gids.items() if len(v) == 1}
        new_gids = set(self.gid_to_node.keys()) - set(d_gid_to_taxonomy.keys())

        for node, taxonomy in self.node_to_taxon_labels.items():
            node_label = self.node_labels[node]

            # Skip those that are not set (i.e. leaf nodes).
            if node_label.bootstrap is None:
                continue

            for taxon in node_label.taxonomy.iter_ranks():
                cur_rank = taxon[0]
                if cur_rank != rank:
                    continue
                cur_rank_parent = RANKS[RANKS.index(cur_rank) - 1]
                cur_rank_child = RANKS[RANKS.index(cur_rank) + 1]

                # Get the median RED dict values for those
                cur_rank_median_red = red_dict.data[cur_rank]
                cur_rank_parent_red = red_dict.data[cur_rank_parent]
                cur_rank_child_red = red_dict.data[cur_rank_child]
                original_node_bs = node_label.bootstrap
                parent_red_thresh = (cur_rank_parent_red + cur_rank_median_red) / 2
                child_red_thresh = (cur_rank_child_red + cur_rank_median_red) / 2
                original_node_red_delta_to_median = abs(cur_rank_median_red - node_label.red)

                # Check up the tree (this is valid for all cases)
                cur_node = node.parent_node
                while cur_node is not None:
                    cur_node_label = self.node_labels[cur_node]
                    cur_node_bs = cur_node_label.bootstrap
                    cur_node_red = cur_node_label.red

                    # Check if we have exceeded the maximum RED range, stop
                    if cur_node_red < parent_red_thresh:
                        break

                    # Check if this node is within the RED range
                    if parent_red_thresh <= cur_node_red <= child_red_thresh:

                        # Check if this node would move the RED value closer to the median
                        cur_node_red_delta_to_median = abs(cur_rank_median_red - cur_node_red)
                        # cur_node_closer_to_median_red = cur_node_red_delta_to_median <= original_node_red_delta_to_median
                        #
                        # # Stop processing it if we are going to bring it further away from the median RED range
                        # if not cur_node_closer_to_median_red:
                        #     break

                        # Now that we have checked the RED criteria, make sure that
                        # this will actually improve the bootstrap support
                        if cur_node_bs >= min_bs or cur_node_bs > original_node_bs:

                            # Verify that we will not be incorporating any nested ranks
                            # by doing this change
                            would_cause_conflicts = False
                            for desc_node in self.node_to_descs[cur_node]:
                                desc_node_tax = self.node_labels[desc_node].taxonomy

                                # This would create a nested rank, we cannot go
                                # any higher
                                if desc_node_tax:
                                    cur_desc_node_tax_for_rank = desc_node_tax.get_rank(rank)
                                    cur_desc_node_tax_for_rank = {x for x in cur_desc_node_tax_for_rank if len(x) > 3}
                                    if len(cur_desc_node_tax_for_rank) > 0 and taxon not in cur_desc_node_tax_for_rank:
                                        would_cause_conflicts = True
                                        break

                            # Stop going higher as this will not improve
                            if would_cause_conflicts:
                                break

                            # Otherwise, we can use this as a candidate node
                            # for the up direction

                            # This would be an improvement to the RED range
                            if cur_node_red_delta_to_median < original_node_red_delta_to_median:
                                out[taxon].append((cur_node, node))

                            # This would move it further from the RED range
                            # we only want to do this if it's going to improve
                            # the bootstrap support
                            else:
                                if node_label.bootstrap is None:
                                    # Don't move it?
                                    print("UNTESTED CODE")
                                else:

                                    # This should only be done if the original
                                    # node was not supported
                                    if node_label.bootstrap < min_bs:
                                        print("UNTESTED CODE")
                                        out[taxon].append((cur_node, node))


                    # Keep going up the tree
                    cur_node = cur_node.parent_node

                # -------------------------------


                # We will not check down the tree for those that taxa that
                # are placed on a leaf node.
                if node_label.bootstrap is None:
                    continue

                # We have a well-supported node
                if node_label.bootstrap >= min_bs:

                    # We do not move well-supported taxa down the tree if it will
                    # cause the RED to move further away from the median
                    if node_label.red > child_red_thresh:
                        continue

                    # We have a well supported node, but moving down would be beneficial
                    else:
                        pass

                # This is not a well-supported node, we should try move it down
                else:
                    pass


                # Get all descendant tax labels (from the metadata)
                d_child_to_tax_labels = defaultdict(set)
                for cur_node_child in node.child_nodes():
                    for cur_node_child_desc in self.node_to_descs[cur_node_child]:
                        if cur_node_child_desc.is_leaf():
                            expected_tax = d_gid_to_taxonomy.get(self.node_labels[cur_node_child_desc].gid)
                            if expected_tax:
                                d_child_to_tax_labels[cur_node_child].update(expected_tax.get_rank(rank))

                would_cause_polyphyly = sum([taxon in x for x in d_child_to_tax_labels.values()]) > 1

                # If we would cause polyphyly by moving the current label down, then we can't do anything
                # although this could be handled in future versions
                if would_cause_polyphyly:
                    continue

                # Since we have made it here, then we should explore what nodes
                # this could be moved to

                # Otherwise, go over the children that wouldn't cause polyphyly
                queue = deque([k for k,v in d_child_to_tax_labels.items() if taxon in v])
                while len(queue) > 0:
                    cur_node = queue.popleft()

                    # As we have made it here, then this is by extension a valid
                    # node to consider, unless it's not well supported or the RED
                    # is terrible
                    cur_node_label = self.node_labels[cur_node]

                    # This is a leaf node, do not go down any further but include it
                    # regardless of the RED range
                    # given that the prior node had bad support
                    if cur_node_label.bootstrap is None:
                        out[taxon].append((cur_node, node))
                        continue

                    else:
                        # If this is well supported then take it
                        if cur_node_label.bootstrap >= min_bs:
                            out[taxon].append((cur_node, node))

                            # Given we have found a good node, we want to stop if
                            # we straying further from the median
                            cur_node_delta_new = abs(cur_rank_median_red - cur_node_label.red)
                            if cur_node_delta_new > original_node_red_delta_to_median:
                                continue

                        # This is an improvement, keep searching
                        elif cur_node_label.bootstrap > node_label.bootstrap:
                            out[taxon].append((cur_node, node))



                    d_node_child_contains_taxon = defaultdict(set)
                    # Go over each child to see what branches (if any) would not cause polyphyly
                    for cur_node_child in cur_node.child_nodes():
                        for cur_node_child_desc in self.node_to_descs[cur_node_child]:
                            if cur_node_child_desc.is_leaf():
                                expected_tax = d_gid_to_taxonomy.get(self.node_labels[cur_node_child_desc].gid)
                                if expected_tax:
                                    d_node_child_contains_taxon[cur_node_child].update( expected_tax.get_rank(rank))

                    # Moving this taxon label down would cause polyphyly
                    would_cause_polyphyly = sum([taxon in x for x in d_node_child_contains_taxon.values()] ) > 1

                    # Not going to handle this case as of now
                    if would_cause_polyphyly:
                        continue

                    # Add the children to the queue that contain the taxon
                    queue.extend([k for k,v in d_node_child_contains_taxon.items() if taxon in v])


        return out





    def process_refinement(self, d_refinment, rank, min_bs, red_dict, metadata):

        if len(d_refinment) == 0:
            return

        # Get the median RED dict values for those
        cur_rank = rank[0]
        # cur_rank_parent = RANKS[RANKS.index(cur_rank) - 1]
        # cur_rank_child = RANKS[RANKS.index(cur_rank) + 1]
        cur_rank_median_red = red_dict.data[cur_rank]
        # cur_rank_parent_red = red_dict.data[cur_rank_parent]
        # cur_rank_child_red = red_dict.data[cur_rank_child]
        # parent_red_thresh = (cur_rank_parent_red + cur_rank_median_red) / 2
        # child_red_thresh = (cur_rank_child_red + cur_rank_median_red) / 2

        # Simply take the nodes that are the most likely
        # 1. Sort by bootstrap >95
        # 2. Sort by red delta to median
        # 3. Sort by bootstrap value desc
        for taxon, node_changes in d_refinment.items():

            changes = dict()
            original_node = [x[1] for x in node_changes][0]
            original_node_red_delta = abs(cur_rank_median_red - self.node_labels[original_node].red)

            for to_node, from_node in node_changes:
                target_node_red_delta = abs(cur_rank_median_red - self.node_labels[to_node].red)

                original_node_bs = self.node_labels[from_node].bootstrap
                target_node_bs = self.node_labels[to_node].bootstrap

                if target_node_bs is None:
                    target_node_bs = 100

                changes[to_node] = (target_node_red_delta, target_node_bs, target_node_bs >= min_bs)
                changes[from_node] = (original_node_red_delta, original_node_bs, original_node_bs >= min_bs)


            changes_sorted = sorted(changes.items(), key=lambda x: (-x[1][2], x[1][0], -x[1][1]))
            best_change_node, _ = changes_sorted[0]

            # Nothing to do here, ideally this wouldn't happen
            assert original_node is not None
            if best_change_node == original_node:
                continue


            # Transfer the taxon label
            logger.info(f'Moving {taxon} from {original_node.label} to {best_change_node.label}')
            transfer_taxon(original_node, best_change_node, taxon)


        return


    def extract_all_taxa_from_tree(self) -> set[str]:
        out = set()
        for node in self.tree.postorder_node_iter():
            x = NodeLabel.from_node(node).taxonomy
            if x:
                for taxon in x.iter_ranks():
                    out.add(taxon)
        return out

