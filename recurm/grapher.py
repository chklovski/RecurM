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
                        return 'Perfect'

        return circularity

    def assign_cluster_label(self, row):
        if row['Circular_Count'] > 0:
            return 'Circular'
        elif row['Imperfect_Count'] == 0:
            return 'Linear'
        else:
            return 'Imperfect'

    def prune_graph(self, graph, min_cluster_size):
        logging.info('Pruning complex graphs.')

        pruned_graphs = []
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


        return pruned_graphs

    def create_graph_from_alignments(self, circular_alignments, linear_alignments, circular_success):
        
        linears = pd.read_csv(linear_alignments, sep='\t', names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'AR', 'LR', 'ANI'])       
        linears['Circularity'] = 'Imperfect'
        
        if circular_success:
            circs = pd.read_csv(circular_alignments, sep='\t', names=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'AR', 'LR', 'ANI'])
            circs['Circularity'] = 'Circular'
            all_alignments = pd.concat([circs, linears])
        else:
            all_alignments = linears

        all_alignments['Circularity'] = all_alignments.apply(lambda row: self.identify_perfect_edge_alignment(row[2], row[7], row[3], row[4], row[8], row[9], row['Circularity']), axis=1)
        all_alignments['Distance'] = all_alignments['LR'] * all_alignments['AR'] * all_alignments['ANI']

        edgeList = all_alignments[[1, 6, 'Distance']].values
        G = nx.Graph()

        for i in range(len(edgeList)):
            G.add_edge(edgeList[i][0], edgeList[i][1], weight=edgeList[i][2])

        return G, all_alignments

    def retrieve_disconnected_subgraphs(self, graph, min_cluster_size, all_alignments):
        initial_subgraphs = (list(graph.subgraph(c) for c in nx.connected_components(graph)))

        sub_graphs = self.prune_graph(initial_subgraphs, min_cluster_size)

        # prune any graphs with nodes < min cluster size
        # also make sure we're not walking through the same sample more than once per cluster
        final_nodes = []
        final_subgraphs = []


        for idx, entry in enumerate(list(sub_graphs)):
            sample_list = []
            for o in entry.nodes():
                sample_list.append(o.split(DefaultValues.DEFAULT_FASTA_HEADER_SEPARATOR)[0])
            # only one contig per sample
            if len(sample_list) == len(set(sample_list)):
                # enough contigs to make valid cluster
                if len(sample_list) >= min_cluster_size:
                    final_nodes.append(entry.nodes())
                    final_subgraphs.append(sub_graphs[idx])

        final_nodes = [list(x) for x in final_nodes]
        if len(final_nodes) == 0:
            logging.warning('No valid clusters found. Check that your --min_contig_length and --min_cluster_size are reasonable.')
            #Probably remove all the intermediate stuff at this point
            sys.exit(0)

        # construct an overview table of results

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


        logging.info('Identified a total of {} clusters. Of these, {} are circular, {} are linear and {} are imperfect.'
                     .format(len(infos),
                             len(infos[infos['Structure'] == 'Circular']),
                             len(infos[infos['Structure'] == 'Linear']),
                             len(infos[infos['Structure'] == 'Imperfect'])))

        return infos, final_subgraphs




