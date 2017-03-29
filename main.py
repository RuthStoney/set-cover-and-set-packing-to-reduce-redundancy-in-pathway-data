
print 'Removing redundancy from pathway sets\n\n'


import meth_functions
functions2= meth_functions


import os
import matplotlib.pyplot as plt
import copy

file_path = os.getcwd()
file_path = file_path.replace('\\','/')
file_path += '/'


user_input = raw_input("to use set packing please type 1 (slow)\nto use standard set cover please press 2\nto use proportional set cover please press 3\nto use hitting set cover please press 4\nto use enrichment set cover please press 5\n")

while user_input not in ['1','2','3','4','5']:
    user_input = raw_input("\n\n sorry I didn't understand your answer\nto use set packing please type 1 (slow)\nto use standard set cover please press 2\nto use proportional set cover please press 3\nto use hitting set cover please press 4\nto use enrichment set cover please press 5\n")

if user_input == '1': set_overlap = 'set_packing'
elif user_input == '2': set_overlap = 'set_cover_unweighted'
elif user_input == '3': set_overlap = 'set_cover_proportional'
elif user_input == '4': set_overlap = 'set_cover_hit_set'
elif user_input == '5': set_overlap = 'set_cover_enriched_paths'
else: print 'broken ' + user_input

if set_overlap == 'set_cover_enriched_paths':
    enrichment_data =1
else:
    enrichment_data =0


  ############################################### OVERLAP #########################################


# FEED IN THE CPDB DATA
user_input = raw_input("do you want to input a pathway textfile? type y for yes\nthe file must be separated by commas or tabs and have .txt extension.\nto use CPDB file press n\n").lower()
while user_input not in ['n','no','yes','y']:
    user_input = raw_input("\n\n sorry I didn't understand your answer\nto you want to input a pathway textfile? type yes\must be separated by commas or tabs and have .txt extension.\nto use CPDB file press n\n")

if user_input == 'y' or user_input == 'yes':
    pathway_genes_dict = functions2.open_pathway_txtfile('text pathways.txt')
else:
    pathway_genes_dict = functions2.CreatePathGeneDict_kegg(file_path+'CPDB_human.tab')
#pathway_genes_dict = {k:pathway_genes_dict[k] for k in pathway_genes_dict.keys()[:200]}

with open('text pathways2.txt', 'w') as f:
    for k,v in pathway_genes_dict.items():
        f.write(k)
        for i in v:
            f.write(',' + i)
        f.write('\n')

mean = functions2.descriptives(pathway_genes_dict, '\nraw pathways')
gene_to_path = functions2.key_val_swapper(pathway_genes_dict)
universe = copy.copy(gene_to_path.keys())

results = functions2.methods_result_scores(pathway_genes_dict)
print 'mean paths per gene raw data: ' + str(results)


if enrichment_data ==1:

    # read in / analyse
    enrich_pvals = functions2.read_enrichment_data()
    enrich_pvals = dict((k, (1.0-v)/1) for (k, v) in enrich_pvals.items())
    pathway_genes_dict = dict((k, v) for (k, v) in pathway_genes_dict.items() if k in enrich_pvals.keys())
    universe = functions2.key_val_swapper(pathway_genes_dict)

    functions2.descriptives(pathway_genes_dict, '\nenriched pathways')
    results = functions2.methods_result_scores(pathway_genes_dict)
    print 'mean paths per gene enriched paths: ' + str(results)

elif set_overlap == 'set_cover_proportion_uncovered':
    enrich_pvals = dict( (k, 1.0/((abs(len(v) - mean)+1)*1000000))    for (k, v) in pathway_genes_dict.items() )



#lens = [len(v) for v in pathway_genes_dict.values()]

# SET COVER
if 'set_cover' in set_overlap: 
    
    overlapping_paths = functions2.overlapping_pathways(pathway_genes_dict)
    gene_to_path = functions2.key_val_swapper(pathway_genes_dict)

    reduced_redundancy_pathways = [] # store the new pathway dictionaries you make

    # get thresholds unless you're doing enriched paths, then it just get complicated with heatmaps
    if set_overlap == 'set_cover_enriched_paths': SCT = [0.0]
    else:SCT = functions2.get_user_thresholds()
    
    # set cover
    for th in SCT: # for each threshold
        ths = str(th)

        # generate the set cover
        if  enrichment_data ==1 or set_overlap == 'set_cover_proportion_uncovered':
            set_cov, cover_paths_order = functions2.set_cover2(gene_to_path, pathway_genes_dict, overlapping_paths, universe, th, set_overlap , enrich_pvals) 
        else:
            set_cov, cover_paths_order = functions2.set_cover2(gene_to_path, pathway_genes_dict, overlapping_paths, universe, th, set_overlap ) 

        # stats on the set cover
        functions2.descriptives(set_cov, '\nset cover' + ths)
        results = functions2.methods_result_scores(set_cov) # gets: abs genes/universe
        print 'mean paths per gene: ' + str(results)

        # generate an out file of the pathways
        if th>0.0: amount_covered = str(int((1-th)*100))
        else: amount_covered = '100'
        functions2.print_dict(set_cov, 'CPDB paths '+ set_overlap + ' ' + amount_covered +' percent genes' + '.txt')
                       
        # store the pathways
        reduced_redundancy_pathways.append([th, set_cov])

        # heat maps for enrichment
        if enrichment_data ==1:
            heat_map_size = 10
            # convert to jpg using http://png2jpg.com/
            functions2.enriched_path_overlap(cover_paths_order, pathway_genes_dict, len(gene_to_path), heat_map_size, 'heatmap after.png')

    if enrichment_data ==1:
        functions2.enriched_path_overlap(enrich_pvals, pathway_genes_dict, len(gene_to_path), heat_map_size, 'heatmap before.jpg')




# SET PACKING
if set_overlap == 'set_packing':
    pathway_genes_dict2 = copy.copy (pathway_genes_dict)
    # write dict of paths that go to R for pathclean_clustR.py
    jaccard_distance_diag_0_bf = functions2.jaccard_simple_matrix (pathway_genes_dict) #jaccard_distance_matrix_bf, jac_matrix_abs_bf,
    gene_list_clean = functions2.genelistmaker(pathway_genes_dict)
 
    set_packing_thresholds = functions2.get_user_thresholds2()
    reduced_redundancy_pathways = functions2.set_packing3(jaccard_distance_diag_0_bf, pathway_genes_dict, len(gene_list_clean), set_packing_thresholds)

    for threshold, set_pack in reduced_redundancy_pathways:
        # stats on the set packing
        functions2.descriptives(set_pack, '\ngenes per path after set packing')
        results = functions2.methods_result_scores(set_pack) # gets: abs genes/universe
        print 'mean paths per gene: ' + str(results)

        # make file of new pathways
        functions2.print_dict(set_pack, 'CPDB paths '+ set_overlap + ' ' + str(threshold) + '.txt')



# make histogram of pathway overlap
jaccard_hist = raw_input("\ndo you want a histogram of jaccard overlap? y or n\n")

while jaccard_hist not in ['y', 'yes', 'n', 'no']:
    jaccard_hist = raw_input("\n\n sorry I didn't understand your answer\ndo you want a histogram of jaccard overlap? y or n\n")


if jaccard_hist == 'y' or jaccard_hist == 'yes':

    plt.close('all')
    if 'set_cover' in set_overlap:
        jaccard_distance_diag_0_bf = functions2.jaccard_simple_matrix(pathway_genes_dict) #, jac_matrix_abs_bf, jaccard_distance_diag_0_bf
        
    jac_before = sorted(list(jaccard_distance_diag_0_bf.values.ravel()))
    #genes_count_paths1 = functions2.gene_path_count(pathway_genes_dict)
    
    for th, pathways in reduced_redundancy_pathways:

        # overlap before and after set packinging
        jaccard_distance_diag_0 = functions2.jaccard_simple_matrix (pathways) #na, na, 
        jac_after_cover = sorted(list(jaccard_distance_diag_0.values.ravel()))
        
        functions2.multi_hist([jac_before, jac_after_cover], 20, ['before redundancy removed', 'after redundancy removed'], 'hist ' + str(th) + '.png', 'log')



