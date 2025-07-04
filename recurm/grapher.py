import networkx as nx
import pandas as pd
import numpy as np
import logging
import sys

from recurm.defaultValues import DefaultValues


class ClusterGraph():
    def __init__(self, threads):
        self.nthreads = threads

    def identify_perfect_edge_alignment(self, alen, blen, astart, aend, bstart, bend, circularity):
        if (blen * 0.99) < alen < (blen * 1.01):  # if within 1% of each other
            if (alen * 0.99) < blen < (alen * 1.01):  # if within 1% of each other
                if astart < (bstart + 1):
                    if aend > (bend - 1):
                        return 'Linear'

        return circularity

    def assign_cluster_label(self, row):
        total_aligns = (row['Imperfect_Count'] + row['Linear_Count']) / 2 #don't care about circulars as they're not relevant
        if total_aligns < 10: #less than 10 contigs
            if total_aligns > 4:
                min_linear_needed = total_aligns - 1
            else:
                min_linear_needed = total_aligns
        else:
            min_linear_needed = 0.91 * total_aligns
        if row['Circular_Count'] > 0 or row['Short_Circular_Count'] > 0:
            return 'Circular'
        elif row['Imperfect_Count'] == 0:
            return 'Linear'
        #allow some leniency in labelling linear clusters - if count is at least 4, and > 90% (or n-1) are linear, it should be linear
        #TODO: ensure the rest are shorter
        elif (row['Linear_Count']/2) > min_linear_needed:
            return 'Linear'
        else:
            return 'Imperfect'

    def prune_graph(self, graph, min_cluster_size):
        logging.info('Pruning complex graphs.')

        pruned_graphs = []

        pruned_count = 0
        for subgraph in graph:

            valid_subgroups = []
            subs = sorted(map(sorted, nx.k_edge_subgraphs(subgraph, k=DefaultValues.EDGE_PRUNE_THRESHOLD)))

            for entry in subs:
                if len(entry) > min_cluster_size:
                    valid_subgroups.append(entry)

            len_dict = {}
            if len(valid_subgroups) == 0:  # either a small cluster or a very weakly linked one.
                pruned_graphs.append(subgraph)
            else:
                if len(valid_subgroups) > 1:
                    pruned_count += 1
                nodes_to_remove = list(subgraph)
                # select the biggest graph resulting from pruning operation
                for group_id, entry in enumerate(valid_subgroups):
                    allcontigs = valid_subgroups[group_id]
                    tempdf = pd.DataFrame(data=allcontigs, columns=['Name'])
                    tempdf['Len'] = tempdf['Name'].apply(lambda x: x.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[1]).astype(int)
                    tempdf.sort_values(by='Len', ascending=False)
                    highestValue = tempdf['Len'].values[0]
                    len_dict[group_id] = highestValue

                max_key = max(len_dict, key=len_dict.get)
                selected_cluster = valid_subgroups[max_key]

                #remove nodes cast aside
                for element in selected_cluster:
                    nodes_to_remove.remove(element)


                unfrozen_subgraph = nx.Graph(subgraph)
                unfrozen_subgraph.remove_nodes_from(nodes_to_remove)
                pruned_graphs.append(unfrozen_subgraph)

        logging.info('Pruned {} graphs.'.format(pruned_count))
        return pruned_graphs

    def create_graph_from_alignments(self, circular_alignments, linear_alignments, circular_success):
        
        linears = pd.read_csv(linear_alignments, sep='\t', names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'AR', 'LR', 'ANI'])       
        linears['Circularity'] = 'Imperfect'
        
        if circular_success:
            circs = pd.read_csv(circular_alignments, sep='\t', names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'AR', 'LR', 'ANI', 'Circularity'])
#            circs['Circularity'] = circs.apply(lambda row:
#                                               'Circular' if row['LR'] >= DefaultValues.FIRST_PASS_AVA_LR_CUTOFF else 'Short_Circular', axis=1)
            all_alignments = pd.concat([circs, linears])
        else:
            all_alignments = linears

        all_alignments['Circularity'] = all_alignments.apply(lambda row: self.identify_perfect_edge_alignment(row[2], row[7], row[3], row[4], row[8], row[9], row['Circularity']), axis=1)
        all_alignments['Distance'] = all_alignments['LR'] * all_alignments['AR'] * all_alignments['ANI']

        #set weight to 0 for short circular alignment in order to not count it as a full alignment for min_contig_count
        all_alignments['Distance'] = all_alignments.apply(lambda row: 0 if row['Circularity'] == 'Short_Circular' else row['Distance'], axis=1)

        edgeList = all_alignments[[1, 6, 'Distance']].values
        G = nx.Graph()

        for i in range(len(edgeList)):
            G.add_edge(edgeList[i][0], edgeList[i][1], weight=edgeList[i][2])

        return G, all_alignments

    def retrieve_disconnected_subgraphs(self, graph, min_cluster_size, all_alignments, outdir, noprune=False):

        initial_subgraphs = (list(graph.subgraph(c) for c in nx.connected_components(graph)))

        if not noprune:
            sub_graphs = self.prune_graph(initial_subgraphs, min_cluster_size)
        else:
            logging.info('Skipping graph pruning due to --noprune flag.')
            sub_graphs = initial_subgraphs

        logging.info(f'Min cluster size: {min_cluster_size}')

        # prune any graphs with nodes < min cluster size
        # also make sure we're not walking through the same sample more than once per cluster
        final_nodes = []
        final_subgraphs = []


        for idx, entry in enumerate(list(sub_graphs)):
            weight_count = 0
            sample_list = []
            for o in entry.nodes():
                sample_list.append(o.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])

            #not counting 0-weighted edges as those are short_circular alignments

            for _, _, data in entry.edges(data=True):
                if data['weight'] > 0:
                    weight_count += 1


            # only one contig per sample
            if len(sample_list) == len(set(sample_list)):
                # enough contigs to make valid cluster
                if len(sample_list) >= min_cluster_size and (weight_count >= min_cluster_size-1):

                    final_nodes.append(entry.nodes())
                    final_subgraphs.append(sub_graphs[idx])

        final_nodes = [list(x) for x in final_nodes]
        if len(final_nodes) == 0:
            logging.warning('No valid clusters found. Check that your --min_contig_length and --min_cluster_size are reasonable.')
            #TODO: Probably remove all the intermediate stuff at this point
            sys.exit(0)

        # TODO: move to after filtering is performed
        # nx.write_gpickle(final_nodes, '{}/{}'.format(outdir, DefaultValues.GRAPH_PICKLE_NAME), protocol=4)

        all_alignments['Key'] = all_alignments[1] + all_alignments[6]
        all_alignments['Key2'] = all_alignments[6] + all_alignments[1]

        circ_dict = dict(zip(all_alignments['Key'], all_alignments['Circularity']))
        circ_dict2 = dict(zip(all_alignments['Key2'], all_alignments['Circularity']))
        circ_dict = {**circ_dict, **circ_dict2}



        list_of_cluster_summaries = []
        id_dict = {}

        for idx, entry in enumerate(final_nodes):
            sizes = []
            samples = []
            members = len(entry)
            circ_count = 0
            short_circ_count = 0
            linear_count = 0
            imprefect_count = 0

            for contig in entry:
                sizes.append(float(contig.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[1]))
                samples.append(contig.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])

            # assess edge quality:
            edges = list(final_subgraphs[idx].edges)
            for e in edges:
                outcome = circ_dict[str(''.join(e))]
                if outcome == 'Linear':
                    linear_count += 1
                elif outcome == 'Circular':
                    circ_count += 1
                elif outcome == 'Short_Circular':
                    short_circ_count += 1
                elif outcome == 'Imperfect':
                    imprefect_count += 1

            l1 = np.array(sizes).mean()
            l2 = np.array(sizes).min()
            l3 = np.array(sizes).max()

            cluster_summary = {'ID': idx,
                 'Average_Bp_Size': l1,
                 'Smallest_Contig_Size': l2,
                 'Largest_Contig_Size': l3,
                 'Total_Contigs': members,
                 'Circular_Count': circ_count,
                 'Short_Circular_Count': short_circ_count,
                 'Linear_Count': linear_count,
                 'Imperfect_Count': imprefect_count}

            df = pd.DataFrame(data=cluster_summary, index=[0])
            list_of_cluster_summaries.append(df)
            id_dict[idx] = entry

        infos = pd.concat(list_of_cluster_summaries)
        infos.reset_index(drop=True, inplace=True)

        #assign label
        infos['Structure'] = infos.apply(lambda row: self.assign_cluster_label(row), axis=1)
        infos['Average_Bp_Size'] = infos['Average_Bp_Size'].astype(int)

        return infos, final_subgraphs




