import re
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pylab
import copy
import time

def open_pathway_txtfile(infile):
    infile = raw_input("whats the name of your file?\n").lower()
    sep = raw_input("press 1 if you file is tab delimited and 2 if your file is comma separated\n").lower()
    while sep not in ['1','2']:
        sep = raw_input("\n\sorry I didnt understand your answer.\npress 1 if you file is tab delimited and 2 if your file is comma separated\n")
    if sep =='1': sep = '\t'
    if sep =='2': sep = ','

    pathways = dict()
    with open(infile) as f:        
        for line in f:
            line = line.split(sep)
            pathways[line[0].strip()] = set()
            for i in line[1:]:
                pathways[line[0].strip()].add(i.strip())
    return pathways


def CreatePathGeneDict_kegg(filename):
    # get tab separated dictionary of pathways and genes
    
    pathway_dict=dict()
    repeated_sets = 0
 #   data_source = dict()
#    data_source_size = dict()


    with open (filename) as f:   #"CPDB_pathways_full_genes .tab"
        for line in f:
            if 'pathway	external_id' in line:
                continue
            
            # get the details from the line
            line = line.split("\t")

            ds = line[2].strip()
##            if database =='ALL' or database == ds:
##                continue
##            if keep_db:
##                if ds not in keep_db: continue

            pathway = re.sub(r'[^\w\s]',' ',line[0]) # remove punctuation
            pathway = pathway.strip()
            pathway = pathway.lower()

            id_list = line[3]
            id_list = id_list.upper()
            id_list = id_list.split(',')
            id_list = [x.strip() for x in id_list]
            id_list = set(id_list)

##            if ds not in data_source.keys():
##                data_source[ds]=1
##                data_source_size[ds]=[len(id_list)]
##                data_source_venn[ds] = set()
##            else:
##                data_source[ds]+=1
##                data_source_size[ds].append(len(id_list))

            if id_list in pathway_dict.values():
                repeated_sets +=1
                continue

            # dont let pathways with identical values in
            if pathway not in pathway_dict.keys():
                pathway_dict[pathway] = set()
            else:
                pathway = pathway + ' ' + ds
                if pathway not in pathway_dict.keys():
                    pathway_dict[pathway] = set()
                else:                   
                    ctr = 2
                    while pathway + ' ' + str(ctr) in pathway_dict.keys():
                        ctr = ctr+1
                    pathway = pathway + ' ' + str(ctr)
                    pathway_dict[pathway] = set()

            for gene in id_list:
                pathway_dict[pathway].add(gene)


    print '\ntotal number of pathways in data set\t' + str(len(pathway_dict))
    print str(repeated_sets) + ' sets skipped because of identical values'
    return pathway_dict


def descriptives(dictionary, title):
    if title: print title
    

    lengthsDict = dict()
    total_for_mean = 0
    allgenes = set()
    for k in dictionary.keys():
        lengthsDict[k]=list()
        lengthsDict[k] = len(dictionary.get(k))
        total_for_mean = total_for_mean + len(dictionary.get(k))
        allgenes = allgenes.union(dictionary.get(k))

    lenlist = sorted(lengthsDict.values())
    if len(lenlist)>0:
        if len(lenlist) % 2 == 0:
            'even case'
            index = (len(lenlist)/2)-1
            median = sum([lenlist[index], lenlist[index+1]])/2.0
        else:
            'odd case'
            median = lenlist[(len(lenlist)/2)]*1.0

        mean =  total_for_mean*1.0/len(lengthsDict.keys())
        
 #       standard_dev(mean, lenlist)
    print 'no. nodes: ' + str(len(dictionary)) + '\ttotal genes/gos ' + str(len(allgenes))
    print 'median: ' + str(median) + ' range: ' + str(lenlist[0]) + ' - ' + str(lenlist[-1] ) + ' mean: ' + str(mean)
    
    return mean


def get_user_thresholds(): # set cover
    user_input = raw_input("\nwould you like the algorithm to stop before all the genes are covered?\ntype n for no or y for yes\n").lower()

    while user_input not in ['n','no','y','yes']:
        user_input = raw_input("\nwould you like the algorithm to stop before all the genes are covered?\ntype n for no or y for yes\n").lower()

    # turn user input as a list of floats
    if user_input == 'y' or user_input == 'yes':
        
        finished = 0
        while finished == 0:
            
            user_input = raw_input("\nplease let us know which thresholds you would like to use\n type the percentage you want, eg type '99' for 99%\nseparate multiple values with a comma, eg '95,99,99.9'\n")
            user_input = "".join(user_input.split()) # remove all whitespace
            user_input = user_input.split(',')
            
            SCT = []
            for i in user_input:
                try:
                    val = float(i)
                except ValueError:
                    val = 'broken'

                if val == 'broken':
                    print 'sorry we didnt understand your answer\n'
                    quit
                else: SCT.append(1.0-(val/100))
                
            if len(SCT) == len(user_input): finished = 1
            
    else: SCT = [0.0]

    return SCT


def get_user_thresholds2():# set packing

    user_input = raw_input("\nwould you like the algorithm to allow some overlap?\ntype n for no or y for yes\n").lower()

    while user_input not in ['n','no','y','yes']:
        user_input = raw_input("\nwould you like the algorithm to allow some overlap?\ntype n for no or y for yes\n").lower()

    # turn user input as a list of floats
    if user_input == 'y' or user_input == 'yes':
        
        finished = 0
        while finished == 0:
            
            user_input = raw_input("\nplease let us know which how much overlap you want to allow (between 0 and 100), \n eg type '10' if pathways can coverlap by up to 10%\nseparate multiple values with a comma, eg 10, 25.5'\n")
            user_input = "".join(user_input.split()) # remove all whitespace
            user_input = user_input.split(',')
            
            SPT = []
            for i in user_input:
                try:
                    val = float(i)
                except ValueError:
                    val = 'broken'

                if val == 'broken':
                    print 'sorry we didnt understand your answer\n'
                    quit
                else: SPT.append(val/100)
                
            if len(SPT) == len(user_input):
                finished = 1
            
    else: SPT = [1.0]

    return SPT


def key_val_swapper(some_dict):
    new_dict=dict()
    for key in some_dict.keys():
        vals = some_dict.get(key)
        for val in vals:
            if val in new_dict.keys():
                new_dict[val].add(key)
            else:
                new_dict[val]={key}
    return new_dict



def methods_result_scores(pathway_dict): 

    abs_genes = [v for vals in pathway_dict.values() for v in vals] 

    return len(abs_genes)*1.0/len(set(abs_genes))



def jaccard_simple_matrix (pathway_dict):           
    valueslist_similarity =[] # used to create a matrix
    valueslist_distance =[]
    valueslist_abs =[]
    
    for path1, genes1 in pathway_dict.items():
        values1 = [len (genes1.intersection(genes2))/float(len(genes1.union(genes2))) for genes2 in pathway_dict.values()]
        valueslist_distance.append(values1)
        values2 = [1-i for i in values1]
        valueslist_similarity.append(values2)
        values3 = [len (genes1.intersection(genes2)) for genes2 in pathway_dict.values()]
        values3.append(-len(genes1))# take away the mach with itself because later it gets summed
        valueslist_abs.append(values3)

    jac_matrix = pd.DataFrame(valueslist_similarity, columns = pathway_dict.keys(), index = pathway_dict.keys())
    jac_matrix_dist = pd.DataFrame(valueslist_distance, columns = pathway_dict.keys(), index = pathway_dict.keys())
    jac_matrix_dist_diag0 = pd.DataFrame(valueslist_distance, columns = pathway_dict.keys(), index = pathway_dict.keys())
    jac_matrix_dist_diag0.values[[np.arange(jac_matrix_dist.shape[1])]*2]=0
    jac_matrix_abs = pd.DataFrame(valueslist_abs, columns = pathway_dict.keys().append('self comp'), index = pathway_dict.keys())
    
    #return jac_matrix_dist, jac_matrix_abs, jac_matrix_dist_diag0
    return jac_matrix_dist_diag0

   
def multi_hist(datasets, bins, labels, outfile, log, title = 'title', lines = 'bars'):
    # histogram with different coloured bars for data sets
    from pylab import xticks
    
    cols = ['red', 'blue', 'limegreen', 'orange', 'dodgerblue', 'crimson']
    if log == 'log': log = True
    else: log = False

    # because in density plots the area under the historgam == 1 rather than the sum of the bars sorting thte y axis is a pain
    # if you set the weight of every data point to (1/no. data_points) then the sum of the bars ==1
    all_weights=[]
    upper_range = 0
    
    for c, d in enumerate(datasets):
        if 0.0 not in d: d.append(0.0)
        w = np.ones_like(d)/float(len(d))
        all_weights.append(w)

        # make sure data is floats
        if type(d[0]) != float:
            datasets[c] = [float(i) for i in d]

        # find the upper range for the bins
        if max(datasets[c])>upper_range: upper_range = max(datasets[c])


    pylab.figure()
    w = np.ones_like(datasets)/float(len(datasets))

    # work out the upper range to 1 significant figures
    n = 1
    upperbin = roundup(upper_range, 2)
    n, bins2, patches = pylab.hist(datasets, bins, weights= all_weights, range=(0.0,upperbin), histtype='bar',  align='mid', color= cols[:len(datasets)], label=labels, normed=False, log=log)

    tickstep = round(upperbin)/bins
    xticks(np.arange(0, upperbin+tickstep, tickstep))
    pylab.legend()
    if log == False: plt.title(title)
    pylab.savefig(outfile)
    plt.close()



def roundup(x, n):
	# rounds a script up to n decimal places
	r = 1
	while x*r<(10**(n-1)):
		r = r*10
	x=x*r
	x = int(math.ceil(x))
	x2 = float('1' + '0'*(len(str(x))-n))
	return (math.ceil(x / x2) * x2)/r


def print_dic_stringval(dic, outfile):
    with open(outfile, 'w') as f:
        for k, v in dic.items():
            f.write(k + '\t' + v + '\n')




def enrichment_data_transcriptsort(my_genes, infile, outfile):

    keepers = []
    genes_found = []
    genes_lost = set()

#   reads in the enriched genes, looks for them in the CPDB pathways
    with open(infile) as f:
        for line in f:
            if 'NA' in line or 'GeneName' in line or 'category' in line: continue
            line = line.split('\t')
            find_gene = line[0].strip()

            if '-' not in find_gene:
                found_gene1 = find_gene_in_list(find_gene, my_genes)
                if found_gene1 == 0:
                    genes_lost.add(find_gene)
                else:
                    found_gene1 = found_gene1.strip()
                    genes_found.append(found_gene1)
                    line[0] = found_gene1
                    keepers.append(line) 

            else:# - means genes in a multi gene complex
                find_genes = find_gene.split('-')
                for find_gene1 in find_genes:
                    find_gene1 =find_gene1.strip()
                    found_gene1 = find_gene_in_list(find_gene1, my_genes)
                    if found_gene1 == 0: genes_lost.add(find_gene1)
                    else:
                        found_gene1 = found_gene1.strip()
                        genes_found.append(found_gene1)
                        line[0] = found_gene1
                        keepers.append(line)

    print 'found ' + str(len(genes_found)*1.0/(len(genes_found)+ len(genes_lost))*100) + '% of genes in OA dataset - ' + str(len(set(genes_found))) + ' unique after matching found'
    
    with open(outfile, 'w') as f:
        for line in keepers:
            for i in line[:-1]:
                f.write(i + '\t')
            f.write(line[-1])
    
    return genes_lost, genes_found



def find_gene_in_list(find_gene, my_genes):
    found = 0
    if find_gene in my_genes:
        found = find_gene
        return found
    
    else:
        # genes may have longer/shorter names than the ones in cpdb
        # change the enrichment name to the cpdb name
        for ref_gene in my_genes:
            if ref_gene[:3] == find_gene[:3]:
                if len(ref_gene)>len(find_gene): long_gene, short_gene = ref_gene, find_gene
                elif len(ref_gene)<len(find_gene): long_gene, short_gene = find_gene, ref_gene
                elif find_gene == ref_gene:
                    found = find_gene
                    return found
                else: continue

                if long_gene[:len(short_gene)] == short_gene:
                    found = ref_gene
                    return found
    return found
                        
   

def read_enrichment_data():

    user_input = raw_input("\ndo you want to input a enrichment textfile? type y for yes\nthe file must be separated by commas or tabs and have .txt extension.\nto use arthritus file press n\n").lower()
    while user_input not in ['n','no','yes','y']:
        user_input = raw_input("\n\n sorry I didn't understand your answer\nto you want to input a pathway textfile? type yes\must be separated by commas or tabs and have .txt extension.\nto use CPDB file press n\n")
        
    if user_input == 'y' or user_input == 'yes':
        enrichment_file = raw_input("whats the name of the enrichment file\n")
        sep = raw_input("press 1 if you file is tab delimited and 2 if your file is comma separated\n").lower()
        while sep not in ['1','2']:
            sep = raw_input("\n\sorry I didnt understand your answer.\npress 1 if you file is tab delimited and 2 if your file is comma separated\n")
        if sep=='1': sep = '\t'
        else: sep = ','
    else:
        enrichment_file = 'enriched p values.txt'
        sep = '\t'
    pval_cutoff = raw_input("what p-value float do you want?\n")
    pval_cutoff = float(pval_cutoff)
   
    
    # reads in files from R goseq expression data
    pvals = dict()

    with open(enrichment_file) as f:
        for line in f:
            path, pval = line.split(sep)
            if float(pval.strip())< pval_cutoff:
                pvals[path.strip()] = float(pval.strip())
    return pvals



def overlapping_pathways(pathway_dict):

    overlap = dict()
    for k1,v1 in pathway_dict.items():
        for k2,v2 in pathway_dict.items():
            if k1 == k2: continue
            
            if len(v1.intersection(v2))>0:
                if k1 not in overlap.keys(): overlap[k1] = set()
                overlap[k1].add(k2)
    return overlap



def gene_path_count(dic):
    d1 = key_val_swapper(dic)
    c = [len(v) for k,v in d1.items()]
    return c





 
def set_cover2(gene_paths, path_genes, overlapping_paths, whole_universe, threshold, ranking_type, extra_weights = []):
    
    if ranking_type == 'set_cover_unweighted': print 'using standard set cover algorithm'
    elif ranking_type == 'set_cover_proportional': print 'using proportional set cover algorithm'
    elif ranking_type == 'set_cover_hit_set': print 'using hitting set cover algorithm'
    elif ranking_type == 'set_cover_enriched_paths': print 'using enriched paths set cover algorithm'
    else: print 'lost'

    universe = set(copy.copy(whole_universe)) # genes left that you need to get
    safe_paths = set()
    safe_order = []

    while len(universe)> len(whole_universe) * threshold:
        
        if len(universe) == len(whole_universe):
            ranked_paths = path_ranking3(gene_paths, path_genes, universe, ranking_type, extra_weights)
        else:
            ranked_paths = path_ranking3(gene_paths, path_genes, universe, ranking_type, extra_weights,  safe_paths, ranked_paths)

        overlap = set()
        last_score = ranked_paths[0][0]
        
        for score, path in sorted(ranked_paths, reverse = True):
            #print path + '\t'+ str(score) + '\t' + str(len(path_genes[path].intersection(universe)))
            # only need to recalculate when the score changes
            if score!= last_score and len(overlap)>0: break 
            last_score = score
            # because paths that have overlap will have a worse score
            if path in overlap: continue

            safe_paths.add(path)
            safe_order.append([score, len(safe_order), path])
            universe = universe - path_genes[path]
            
            if len(universe) ==0: quit

            # if the safe path has any others that overlap get them
            if path in overlapping_paths.keys():
                overlap = overlap.union(overlapping_paths[path])
                overlap = overlap - safe_paths

    safe_dict=dict()
    for i in safe_paths:
        safe_dict[i]= path_genes[i]
                    
    return safe_dict, safe_order #, overlap_scatter
        

def path_ranking3(gene_paths, path_genes, universe, score_type, pvals_med = None,  saved_paths=[], ranked_paths=[]):
    # # (needed_genes / not_needed_genes) / square root(abs(path_length - median) +1)
 
    # delete the saved paths and the paths that need recalculating
    if len(saved_paths)>0:
        # print 'paths going in \t' + str(len(dev_paths1))
        calc_paths = set()
        for score, path in ranked_paths:
            if path in saved_paths: continue
            else: calc_paths.add(path)
    else:
        calc_paths = path_genes.keys()
 
    new_ranking = list()
    for path in calc_paths:
 
        score1 = len(path_genes[path].intersection(universe))
                  
        if score_type == 'set_cover_unweighted':
            new_ranking.append([score1, path])

        elif score_type == 'set_cover_enriched_paths':
            if score1>0: score1 = 1.0
            else: score1 = 0
            new_ranking.append([score1 * pvals_med[path], path])
 
        elif score_type == 'set_cover_proportional':
            unweighted_score = score1*1.0/len(path_genes[path])
            if pvals_med: unweighted_score = unweighted_score + pvals_med[path]
            new_ranking.append([unweighted_score, path])
 
        elif score_type == 'set_cover_hit_set':
            score2 = [1.0/len(gene_paths[gene]) for gene in path_genes[path] if gene in universe]
            score2 = sum(score2)/len(path_genes[path])
            new_ranking.append([score2, path])

##        elif score_type == 'set_cover_unique_genes':
##            score2 = [gene for gene in path_genes[path] if len(gene_paths[gene]) == 1]
##            if score2:
##                score2 = len(score2)*1.0/len(path_genes[path])
##            else:
##                score2 = 0
##            new_ranking.append([score1+score2, path])


        else: print  'ranking didnt happen'

    return sorted(new_ranking, reverse=True)
 


def enriched_path_overlap(scores, genes, total_genes, heat_map_size, outfile):
    sys.path.append('/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
    import image
    
    if type(scores) == dict:
        scores = sorted([[v,v,k] for k,v in scores.items()], reverse=True)
        
    hm = np.ndarray(shape=(heat_map_size,heat_map_size))
    labels = []
    saved_genes = set()

    for c1, path1 in enumerate(scores[:heat_map_size]):
        labels.append(path1[2].replace('homo sapiens  human', '').strip())
        saved_genes = saved_genes.union(genes[path1[2]])
        for c2, path2 in enumerate(scores[:heat_map_size]):
            i = len(genes[path1[2]].intersection(genes[path2[2]]))*1.0/len(genes[path1[2]])
            if c1==c2: hm[c1, c2] = None
            else:
                hm[c1, c2] = i
                #if i == 1.0: print path1[2] + ' is a subset of ' + path2[2]

#   do plot
    fig = plt.figure(figsize=(8, 6)) 
    ax = fig.add_subplot(1,1,1)
    cax = ax.imshow(hm,  interpolation='none' , cmap='coolwarm', vmin=0, vmax=1)

    # marks
    for x in range(0,heat_map_size):
        for y in range(0,heat_map_size):
            score = round(hm[y,x] *100)/100
            if score>0.0:
                plt.annotate( score , xy = (x-0.15,y), size = 9)
    
    
    # sort labels
    plt.yticks(np.arange(0,heat_map_size, 1))
    plt.xticks(np.arange(0,heat_map_size, 1))
 #   labels = [' '*(45-len(i)) + i for i in labels]
    ax.set_xticklabels(np.asarray(labels), rotation=25, horizontalalignment = 'right', fontsize=9)
    ax.set_yticklabels(np.asarray(labels), fontsize=12)

    # save and close
    fig.colorbar(cax, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
    
    plt.savefig(outfile, bbox_inches='tight')
    image.open('testplot.png')#.save(outfile.replace('.png', '.jpg',),'JPEG')
    plt.close('all')

    print str(len(saved_genes)*1.0/total_genes) + ' genes covered by top '+ str(heat_map_size) +' pathways'
           



def genelistmaker(finalDict):
    vals=set()
    for k in finalDict.keys():
        v = finalDict.get(k)
        for i in v:
            if i.upper().strip() not in vals:
                vals.add(i)
    print str(len(vals)) + ' unique values'
    return vals



def set_packing3(jaccard_distance_matrix, path_genes, max_genes, ks): # Independent set (graph theory) path_genes, max_genes
    print '\nstarted set packing at ' + time.strftime('%X %x')
    
    # rank pathways in order of their % overlap with other pathways
    jac_distance_matrix = copy.deepcopy(jaccard_distance_matrix)
    paths = jac_distance_matrix.index.values
    start_no = len(paths)
    return_pathways = []
 #   scorelist=[]
 
    for k in ks:
        jac_distance_matrix = copy.deepcopy(jaccard_distance_matrix)
        keepers = dict()
        lost_paths = set()
        safe_genes = set()

        while len(keepers) + len(lost_paths) < start_no:
            #print start_no - len(keepers) + len(lost_paths)
            counter=0
            recalc = 0
            overlap = rank_paths(paths, jac_distance_matrix)

            while counter< len(overlap) and recalc == 0:
            #for score, path in overlap:

                score, path = overlap[counter]

                # skip if the path is in lost paths or already in keepers
                if path in lost_paths or path in keepers.keys():
                    counter+=1
                    continue
                if not path: print str(len(overlap)) + '\t' + str(counter) + '\t' + ' path disapeared'

                # keep the path and get its genes
                keepers[path] = path_genes.get(path)
                try:
                    safe_genes = safe_genes.union(path_genes.get(path))
                except Exception as e:
                    if path not in path_genes.keys(): print(path + 'not in dict')

                # find the overlapping paths
                lost_neighbours = []
                if score > 0.0:
                    overlap_paths = jac_distance_matrix.loc[path]
                    lost_neighbours = [overlap_paths.index[c] for c, i in enumerate(overlap_paths) if i>k]
                    
                    for n1 in lost_neighbours:
                        if n1 in keepers: print 'super weird jaccard is symetrical'
                        else: lost_paths.add(n1)

                counter+=1
                
                # remove the overlap paths from the distance matrix
                if lost_neighbours:
                    for overlap_p in lost_neighbours:
                        jac_distance_matrix = jac_distance_matrix.drop(overlap_p, 0)
                    for overlap_p in lost_neighbours:
                        jac_distance_matrix = jac_distance_matrix.drop(overlap_p, 1)
                    paths = jac_distance_matrix.index.values

                    #rerank the scores
                    recalc = 1
                else: continue

        print '\n'+str(len(safe_genes)*1.0/max_genes) + ' genes saved by set packing, '  + str(len(keepers)*1.0/start_no) +' paths kept at k of '+ str(k)
        descriptives(keepers, '\nset packing ' + str(k))
        gene_to_path = key_val_swapper(keepers)
 
        results = methods_result_scores(keepers)
        print ('mean number of pathways per gene: ' + str(results))
        return_pathways.append([k, keepers])


    print ' finished set packing at ' + time.strftime('%X %x') + '\n'
    return return_pathways




def rank_paths(paths, jaccard_distance_matrix):
    # textbook weight(1)/neighbourhood size
    overlap = []

    for path in paths:
        col = jaccard_distance_matrix.loc[path]
        #overlap.append([sum(col), path])
        col = [i for i in col if i > 0.0]
        overlap.append([len(col), path])

    return sorted(overlap)

          
def dict_length(dic):
    dic_length = dict()

    for k,v in dic.items():
        dic_length[k]=len(v)
    return dic_length
        

def print_dict(dic, outfile):
    with open(outfile, 'w') as f:
        for k, v in dic.items():
            f.write(k)
            for i in v:
                f.write('\t' + i)
            f.write('\n')



def get_hist_filename():
    # get a file name
    outfile = raw_input("\nplease provide a file name ending in .png or .pdf or .jpg\n")
    filename, filetype = outfile.split('.')
    while not filename:
        jaccard_hist = raw_input("\n\n sorry I didn't understand your answer1\nplease provide a file name ending in .png or .pdf or .jpg\n")
        filename, filetype = outfile.split('.')

    while filetype not in ['png', 'pdf', 'jpg']:
        jaccard_hist = raw_input("\n\n sorry I didn't understand your answer2\nplease provide a file name ending in .png or .pdf or .jpg\n")
        filename, filetype = outfile.split('.')
    return outfile


