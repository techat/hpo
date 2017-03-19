import pandas as pd
import pprint

HPO_FILE = 'all_hpo_terms_and_synonyms.txt'
PHENOTYPE_FILE = 'sample_patient_phenotype.txt'
PHENOTYPE_TO_GENE_FILE = 'Expanded_OMIM_ALL_FREQUENCIES_phenotype_to_genes_without_synonym.txt'
#PHENOTYPE_TO_GENE_FILE = 'Expanded_OMIM_FREQUENT_FEATURES_phenotype_to_genes_without_synonym.txt'

# In practice, we don't need to screen through all genes. Only a few gene candidates are examined.
#CANDIDATE_GENES = ['PTPN11', 'BRCA1']

CANDIDATE_GENES = ['ANK1']
#CANDIDATE_GENES = ['GRIN1', 'HDGFRP2']
#CANDIDATE_GENES = ['GNAS']

## Get levels of each hpo term; e.g., HP:0000001 is the first level; its direct children are the second level.
## Phenotypic abnormality is level 1; Abnormality of prenatal development or birth is level 2.
## In the search for the common ancestors, we should stop at level 2; if two terms have the common ancestor below or equal at level 2, then we deem them to be semantically similar.

hpo_subclass = dict()
with open('hpo_subclasses.txt', 'rb') as f:
    # skip the first header line
    f.readline()
    for line in f.readlines():
        line = line.rstrip()
        parts = line.split('\t')
        hpoid = parts[0]
        subclasses = parts[1].strip('[]').replace("'", "")
        subclasses = [_.strip() for _ in subclasses.split(',')]
        hpo_subclass[hpoid] = subclasses

hpo_levels = dict()
level = 0
root = 'HP:0000001'
curr_level_nodes = [root] 
while curr_level_nodes:
    next_level_nodes = []
    for node in curr_level_nodes:
        hpo_levels[node] = level
        if node in hpo_subclass:
            next_level_nodes += hpo_subclass[node]
    curr_level_nodes = next_level_nodes
    level += 1
     

hpo_superclass = dict()
with open('hpo_superclasses.txt', 'rb') as f:
    # skip the first header line
    f.readline()
    for line in f.readlines():
        line = line.rstrip()
        parts = line.split('\t')
        hpoid = parts[0]
        superclasses = parts[1].strip('[]').replace("'", "")
        superclasses = [_.strip() for _ in superclasses.split(',')]
        hpo_superclass[hpoid] = superclasses

hpo_name = dict()
with open('hpo_names.txt', 'rb') as f:
    # Skip the first header line
    f.readline()
    for line in f.readlines():
        line = line.rstrip()
        parts = line.split('\t')
        hpoid = parts[0]
        hponame = parts[1]
        hpo_name[hpoid] = hponame        

from collections import deque
tree = hpo_superclass
def getAllAncestorsBFS(node, tree):
    P, Q = [node], deque([node])
    while Q:
        u = Q.popleft()
        if u not in tree:
            continue
        for v in tree[u]:
            if v in P: continue
            if v not in tree: continue
            P.append(v)
            Q.append(v)
    # Append All ('HP:0000001') as the highest ancestor 
    P.append('HP:0000001')
    return P 

def findLowestCommonAncestor(node1, node2, tree):
    node1_ancestors = getAllAncestorsBFS(node1, tree) 
    node2_ancestors = getAllAncestorsBFS(node2, tree)
    common_ancestors = [] 
    for ancestor in node2_ancestors:
        if ancestor in node1_ancestors:
            name = hpo_name[ancestor]
            level = hpo_levels[ancestor]
            common_ancestors.append((ancestor, level, name)) 
    lca = common_ancestors[0]
    return lca, common_ancestors

## Generate antonyms of a word using NLTK and wordnet
from nltk.corpus import wordnet as wn
def get_antonyms(word):
    antonyms = []
    for synset in wn.synsets(word):
        for lemma in synset.lemmas():
            #print synset, lemma, lemma.name(), lemma.antonyms(), synset.similar_tos()
            antonyms += lemma.antonyms()
            similar_synsets = synset.similar_tos() 
            for similar_synset in similar_synsets:
                for lemma_ in similar_synset.lemmas():
                    #print (similar_synset, lemma_, lemma_.name(), 
                    #      lemma_.antonyms(), similar_synset.similar_tos())
                    antonyms += lemma_.antonyms()
    direct_antonyms = list(set(antonyms))
    # Get similar words of the direct antonyms
    for antonym in direct_antonyms:
        #print antonym
        similar_synsets = antonym.synset().similar_tos()
        for synset in similar_synsets:
            for lemma in synset.lemmas():
                #print (similar_synset, lemma, lemma.name(), 
                #      lemma.antonyms(), similar_synset.similar_tos())
                antonyms.append(lemma)
    antonyms = list(set(antonyms))
    antonyms = set([_.name() for _ in antonyms])
    return antonyms

## Map patient phenotypes to hpo standard terms 
with open(PHENOTYPE_FILE, 'rb') as infile:
    phenos = []
    for line in infile:
        if line.startswith('#') or not line.strip():
            continue
        line = line.rstrip()
        phenos += line.split(',')
    phenos = [_.strip() for _ in phenos]

#df_hpo = pd.read_csv(HPO_FILE, sep = '\t')
###df_hpo['hpoid'] = df_hpo['hpoid'].str.split('-').str.get(0)
#hpo_name2id = df_hpo.set_index('hponame').to_dict()['hpoid']
#hpo_terms = hpo_name2id.keys()

with open(HPO_FILE, 'rb') as infile:
    hpo_name2id = dict()
    hpo_id2synonyms = dict()
    # skip the first line
    infile.readline()
    for line in infile.readlines():
        line = line.rstrip()
        parts = line.split('\t')
        hpoid, hponame = parts[0], parts[1]
        hpo_name2id[hponame] = hpoid
        hpoid = hpoid.split('-')[0]
        if hponame in hpo_id2synonyms:
            hpo_id2synonyms[hpoid].append(hponame)
        else:
            hpo_id2synonyms[hpoid] = [hponame]
    hpo_terms = hpo_name2id.keys()  


import re
from nltk.stem.wordnet import WordNetLemmatizer
# Use re for split hpo terms
# Use lemmatizer to lemmatize words e.g., disturbances --> distrubance

def isOppositeMeaning(words_in_pheno_only, words_in_hpo_only):
    # If any word ends with ness, convert to its noun form
    words_in_pheno_only |= set([_[0:-4] for _ in words_in_pheno_only if _.endswith('ness')])
    words_in_hpo_only |= set([_[0:-4] for _ in words_in_hpo_only if _.endswith('ness')])
    # corner case: gain
    if 'gain' in words_in_hpo_only:
        words_in_hpo_only.add('increased')
        words_in_hpo_only.add('steady')
    if 'gain' in words_in_pheno_only:
        words_in_pheno_only.add('increased')
        words_in_pheno_only.add('steady')
    # Poor weight gain and weight gain are in opposite meaning
    if (words_in_pheno_only == set(['poor']) and words_in_hpo_only == set()) or (words_in_hpo_only == set(['poor']) and words_in_pheno_only == set()):
        print "Opposite meaning due to word 'poor'", words_in_pheno_only, words_in_hpo_only
        return True
    for word in words_in_pheno_only:
	try:
	    antonyms = get_antonyms(word)
	except UnicodeDecodeError as e:
	    continue
	if antonyms.intersection(words_in_hpo_only):
            print "Opposite meaning: ", antonyms.intersection(words_in_hpo_only)
	    return True  
    for word in words_in_hpo_only:
	try:
	    antonyms = get_antonyms(word)
	except UnicodeDecodeError as e:
	    continue
	if antonyms.intersection(words_in_pheno_only):
            print "Opposite meaning: ", antonyms.intersection(words_in_hpo_only)
	    return True  
    return False 


# If the opposite meanings are found for a hpoid, then its synonyms are all in opposite meanings
unmatches_due_to_opposite_meanings = []

def map2hpo(pheno):
    global unmatches_due_to_opposite_meanings
    matches = []
    # Use lemmatizer to lemmatize words e.g., disturbances --> distrubance
    lemmatizer = WordNetLemmatizer()
    for hpo in hpo_terms:
        hpoid = hpo_name2id[hpo]
        pheno_ = pheno.lower() 
        hpo_ = hpo.lower()
        if hpo_.startswith('obsolete '):
            continue 
        if pheno_ == hpo_:
            matches.append((hpoid, hpo, 'EXACT', 1.0)) 
            continue
        pheno_words = [_.strip(",()'") for _ in re.split(' |/', pheno_)]
        hpo_words = [_.strip(",()'") for _ in re.split(' |/', hpo_)]

        #print hpo_words
        try:
            pheno_words = [lemmatizer.lemmatize(word) for word in pheno_words]
            hpo_words = [lemmatizer.lemmatize(word) for word in hpo_words]
        except UnicodeDecodeError as e:
            pass
        if pheno_words == hpo_words:
            matches.append((hpoid, hpo, 'LEMMA_MATCH', 1.0)) 
            continue
        if set(pheno_words) == set(hpo_words):
            matches.append((hpoid, hpo, 'REORDER', 1.0)) 
            continue

        common_words = set(pheno_words) & set(hpo_words)
        #if len(common_words) > 0:
        #    print pheno_words, hpo_words
        words_in_pheno_only = set(pheno_words) - common_words
        words_in_hpo_only = set(hpo_words) - common_words
        # check if the pheno and hpo are in opposite meaning
        # if any word in pheno is an antonym of a word in hpo, or the other way
        # we skip this hpo
        sim = float(len(common_words)) / max(len(hpo_words), len(pheno_words))  
        if (sim >= 0.3333) and (pheno_ in hpo_) and (len(hpo_words) <= 10 * len(pheno_words)):
            if not isOppositeMeaning(words_in_pheno_only, words_in_hpo_only):
                matches.append((hpoid, hpo, 'SUB', sim))
            else:
                unmatches_due_to_opposite_meanings.append(hpoid.split('-')[0])
            continue
        if (sim >= 0.3333) and (hpo_ in pheno_) and (len(pheno_words) <= 10 * len(hpo_words)): 
            if not isOppositeMeaning(words_in_pheno_only, words_in_hpo_only):
                matches.append((hpoid, hpo, 'SUPER', sim))
            else:
                unmatches_due_to_opposite_meanings.append(hpoid.split('-')[0])
            continue
        if (sim >= 0.5) and len(pheno_words) >= 2:
            if not isOppositeMeaning(words_in_pheno_only, words_in_hpo_only):
                matches.append((hpoid, hpo, 'OVERLAP', sim))
            else:
                unmatches_due_to_opposite_meanings.append(hpoid.split('-')[0])
    matches_res = []
    for match in matches:
        if match[0].split('-')[0] not in unmatches_due_to_opposite_meanings:
            matches_res.append(match)
    matches = matches_res
    return matches

def filterMatchesOnCommonAncestors(matches):
    ## e.g., for "Central hypotonia", matches is like: 
    ## [('Central', 'SUPER', 0.5), ('Central hypotonia', 'EXACT', '1.0'), ('Hypotonia', 'SUPER', 0.5)]
    ## lca is like ('HP:0005872', 9, 'Brachytelomesophalangy')
    ## The function returns like: [('Stereotypical hand wringing', 'HP:0012171', 0.6666666666666666)]
    if not matches:
        return matches 
    best_match_hpo = sorted(matches, key = lambda x: x[3], reverse = True)[0][0]
    final_matches = []
    for match in matches:
        term = match[0]
        # Remove 'synonym' from the hpoid (term) e.g., HP:0004414-synonym --> HP:0004414
        term = term.split('-')[0]
        lca, common_ancestors = findLowestCommonAncestor(term, best_match_hpo.split('-')[0], hpo_superclass)
        if lca[1] >= 2:
            final_matches.append((match[1], match[0], match[3])) 
    return final_matches    

def map2hpoWithPhenoSynonyms(pheno):
    global unmatches_due_to_opposite_meanings
    matches = map2hpo(pheno)
    direct_matches = filterMatchesOnCommonAncestors(matches)  
    final_matches = direct_matches 
    for match in direct_matches:
        hponame = match[0]
        hpoid = match[1].split('-')[0]
        sim = match[2]
        if sim >= 0.75:
            hposynonyms = hpo_id2synonyms[hpoid]
            for synonym in hposynonyms:
                if synonym != hponame:
                    print synonym
                    matches_ = map2hpo(synonym)
                    indirect_matches = filterMatchesOnCommonAncestors(matches_)
                    final_matches = final_matches + indirect_matches
    final_matches = list(set(final_matches))
    matches_res = []
    for match in final_matches:
        if match[0].split('-')[0] not in unmatches_due_to_opposite_meanings:
            matches_res.append(match)
    final_matches = matches_res
    unmatches_due_to_opposite_meanings = []
    return final_matches

hpoid2gene = dict()
with open(PHENOTYPE_TO_GENE_FILE, 'rb') as f:
    # skip the first header line
    f.readline()
    for line in f.readlines():
	line = line.rstrip()
	parts = line.split('\t')
	geneid, genesymbol, hponame, hpoid = parts[0], parts[1], parts[2], parts[3]
	if hpoid in hpoid2gene:
	    hpoid2gene[hpoid].append(genesymbol)   
	else:
	    hpoid2gene[hpoid] = [genesymbol]

def map2gene(final_matches):
    # Map to genes for each phenotype (one phenotype at one time)
    # Pre-process the final_matches by remove '-synonym' from the hpoid and then unique the hpoids
    final_matches = [_[1] for _ in final_matches] 
    final_matches = [_.split('-')[0] for _ in final_matches] 
    final_matches = list(set(final_matches))
    mapped_genes = []
    for hpoid in final_matches:
        mapped_genes += hpoid2gene[hpoid] 
    ## For each phenotype in patient record, we generated multiple phenotype keywords 
    ## to maximize the mapping to the genes. However, we only count it once for each 
    ## phenotype even if the multiple keywords lead to multiple mappings.
    ## In practice, we don't need to screen through all genes. Only a few gene candidates are examined.
    mapped_genes = list(set(mapped_genes) & set(CANDIDATE_GENES))
    return mapped_genes

def map2geneWithSim(final_matches):
    # Map to genes for each phenotype (one phenotype at one time) with similarity information (word len)
    # final_matches like: [('Postnatal macrocephaly','HP:0005490',0.5), ('Macrocephaly,relative','HP:0004482-synonym', 0.5), ...]
    # map2gene function returns a list while this function returns a dict with genesymbol as key and similarity as value 
    hpoid_sim = dict() 
    for match in final_matches:
        hpoid = match[1].split('-')[0]
        sim = match[2]
        if hpoid not in hpoid_sim or sim > hpoid_sim[hpoid]:
            hpoid_sim[hpoid] = sim
    mapped_genes_score = dict()
    for hpoid in hpoid_sim:
        mapped_genes = hpoid2gene[hpoid]
        sim = hpoid_sim[hpoid]
        for gene in mapped_genes:
            # In practice, we don't need to screen through all genes. 
            # Only a few gene candidates are examined.
            if gene not in CANDIDATE_GENES:
                continue
            if gene not in mapped_genes_score or sim > mapped_genes_score[gene]:
                mapped_genes_score[gene] = sim
    return mapped_genes_score 

from collections import Counter
all_mapped_genes = []
all_mapped_genes_score = {}
gene_phenos = {} 
for pheno in phenos:
    print ("===========================================================================================")
    pprint.pprint(pheno)
    matches =  map2hpo(pheno) 
    #pprint.pprint(matches)
    #final_matches = filterMatchesOnCommonAncestors(matches)
    final_matches = map2hpoWithPhenoSynonyms(pheno)
    pprint.pprint(final_matches)
    mapped_genes = map2gene(final_matches)
    pprint.pprint(mapped_genes)
    all_mapped_genes += mapped_genes
    mapped_genes_score = map2geneWithSim(final_matches)
    for gene in mapped_genes_score:
        if gene not in all_mapped_genes_score:
            all_mapped_genes_score[gene] = mapped_genes_score[gene]
        else:
            all_mapped_genes_score[gene] += mapped_genes_score[gene]
        if gene not in gene_phenos:
            gene_phenos[gene] = [pheno]
        else:
            gene_phenos[gene].append(pheno)

all_mapped_genes_score = [(k,v) for k, v in all_mapped_genes_score.items()]   
all_mapped_genes_score = sorted(all_mapped_genes_score, key = lambda pair: pair[1], reverse = True)

count = Counter(all_mapped_genes)
count = sorted(count.items(), key = lambda pair: pair[1], reverse = True)
print "********************************************************************************************"
print "********************************************************************************************"
print "Count of matches to phenos: "
pprint.pprint(count[0:10])
#print len(count)
print "********************************************************************************************"
print "********************************************************************************************"
print "Score of matches to phenos: "
pprint.pprint(all_mapped_genes_score[0:10])
#print len(all_mapped_genes_score)
print "********************************************************************************************"
print "********************************************************************************************"
print "Matched phenos: "
pprint.pprint(gene_phenos)

'''
Use spacy or other POS tagger to recognize the important NOUN and then do the match; this may help to remove the wrong matches like "Central hypotonia" to "Central"  
Macrocephaly
[('Congenital macrocephaly', 'SUB', 0.5), ('Macrocephaly', 'EXACT', '1.0'), ('Macrocephaly at birth', 'SUB', 0.3333333333333333), ('Macrocephaly, postnatal', 'SUB', 0.5), ('Macrocephaly, progressive', 'SUB', 0.5), ('Macrocephaly, relative', 'SUB', 0.5), ('Postnatal macrocephaly', 'SUB', 0.5), ('Progressive macrocephaly', 'SUB', 0.5), ('Relative macrocephaly', 'SUB', 0.5)]
Central hypotonia
[('Central', 'SUPER', 0.5), ('Central hypotonia', 'EXACT', '1.0'), ('Hypotonia', 'SUPER', 0.5)]
Ventriculomegaly
[('Mild fetal ventriculomegaly', 'SUB', 0.3333333333333333), ('Progressive ventriculomegaly', 'SUB', 0.5), ('Ventriculomegaly', 'EXACT', '1.0')]
Mild global developmental delay
[('Developmental delay', 'SUPER', 0.5), ('Global developmental delay', 'SUPER', 0.75), ('Mild global developmental delay', 'EXACT', '1.0')]
'''

