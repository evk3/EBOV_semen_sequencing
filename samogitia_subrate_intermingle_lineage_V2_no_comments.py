import argparse
import re
import datetime as dt
import baltic as bt
import sys
import collections
import numpy as np
from collections import OrderedDict

#In v2.0, incorporates regex code from S. Whitmer and started reporting info about whether subgroups are monophyletic
#In lineage_trace version, working to incorporate unique branch info from Gytis.

#Intermingle_lineage_V1 is designed to seperate shared lineages into blood-specific and semen-specific 
#lineages, with no shared branches.  Problem with this verion, because it cannot calculate an acute-only
#background rate.  The reason for this is because I did not modify the acute only rate estimate and I did not
#make any efforts to include the eye rate.

#Intermingle_lineage_V2 is designed to also calculate the rate for the branch leading to the 
#eye sample.  Hard-coded this, because a little bit faster than modifying the arguments.  Also,
#fixed the acute-only rate to be all_rates - persistent_rates. yeay!

def overlap(a,b):
    """
    Return the elements shared by two lists in the following format:
    [overlap],[list 1 remainder],[list 2 remainder]
    """
    a_multiset = collections.Counter(a)
    b_multiset = collections.Counter(b)

    overlap = list((a_multiset & b_multiset).elements())
    a_remainder = list((a_multiset - b_multiset).elements())
    b_remainder = list((b_multiset - a_multiset).elements())

    return overlap, a_remainder, b_remainder

def subclade_branch_info(node):
    lengths=0
    numsubs=0
    for subnode in node.children:
        #Needed for starting at the root, something to do with node 0 not having a rate trait
        if 'rate' in subnode.traits:
            lengths+=subnode.length
            numsubs+=subnode.length*subnode.traits['rate']
        if isinstance(subnode,bt.node):
            sublen, subsubs = subclade_branch_info(subnode)
            lengths+=sublen
            numsubs+=subsubs
    return lengths, numsubs

def lineage_trace_info(tMRCA_node, leaves):
    #This function should trace the branches from the leaf to the parent node.
    #The loop should stop tracing when we reach tMRCA node.
    #The traced path is stored in lineage_trace[].
    #We do not want to store any (duplicate) paths traced.
    ll.renameTips(tips)
    for k in ll.Objects: #Iterate over branches
        if isinstance(k, bt.leaf):
            for s in leaves:
                if s == k.numName:
                    s = k  #Assign the branch to the interator.
                    #print "Entering the lineage trace."
                    #print "Initial iterator value is: %s" %(s)
                    while s != tMRCA_node:
                        #print "Lineage trace from: %s" %(k)
                        #print "Current branch is: %s" %(s)
                        if s not in lineage_trace:
                            lineage_trace.append(s)
                            #print "Added %s to lineage trace" %(s)
                        else:
                            #print "Not added to lineage, duplicate branch: %s" %(s)
                            #print k.numName
                            pass
                        s = s.parent
                        #print "Current branch iterates up to: %s" %(s)
                        #print "tMRCA node is: %s" %(tMRCA_node)
        
    return set(lineage_trace)

#Increase recursion limit before starting script
sys.setrecursionlimit(50000)


############ arguments
austechia = argparse.ArgumentParser(description="samogitia.py analyses trees drawn from the posterior distribution by BEAST.\n")

austechia.add_argument('-b','--burnin', default=0, type=int, help="Number of states to remove as burnin (default 0).\n")
austechia.add_argument('-c','--calibrate', default=False, type=bool, help="Set the tree to have absolute time information (default True). Should only be used if tip names contain information about *when* each sequence was collected.\n")
austechia.add_argument('-t','--treefile', type=open, help="File with trees sampled from the posterior distribution (usually with suffix .trees).\n")
austechia.add_argument('-a','--analyses', type=str, nargs='+', help="Analysis to be performed, can be a list separated by spaces.\n")
austechia.add_argument('-o','--output', type=argparse.FileType('w'), default='samogitia.out.txt', help="Output file name (default samogitia.out.txt).\n")
austechia.add_argument('-s','--states', type=str, default='0-inf', help="Define range of states for analysis.\n")
austechia.add_argument('--substr', type=str, nargs='+', help="List separated by spaces. Strings to use to identify names of leafs in subrate groups\n")

args = vars(austechia.parse_args())
burnin, treefile, analyses, outfile, calibration, states, substr = args['burnin'], args['treefile'], args['analyses'], args['output'], args['calibrate'], args['states'], args['substr']

lower,upper=states.split('-')
lower=int(lower)
if upper=='inf':
    upper=np.inf
else:
    upper=int(upper)

try:
    for line in open('banner_samogitia.txt','r'):
        sys.stderr.write('%s'%(line))
except:
    pass

## keeps track of things in the tree file at the beginning
plate=True 
taxonlist=False
treecount=0 ## counts tree
## 

tips={} ## remembers tip encodings

############################# progress bar stuff
Ntrees=10000 ## assume 10 000 trees in posterior sample
barLength=30
progress_update=Ntrees/barLength ## update progress bar every time a tick has to be added
threshold=progress_update ## threshold at which to update progress bar
processingRate=[] ## remember how quickly script processes trees
#############################

available_analyses=['treeLength','RC','Sharp','tmrcas','transitions','subtrees', 'subrate'] ## analysis names that are possible

assert analyses,'No analyses were selected.'
for queued_analysis in analyses: ## for each queued analysis check if austechia can do anything about them (i.e. whether they're known analysis types)
    assert queued_analysis in available_analyses,'%s is not a known analysis type\n\nAvailable analysis types are: \n* %s\n'%(queued_analysis,'\n* '.join(available_analyses))

begin=dt.datetime.now() ## start timer

for line in treefile: ## iterate through each line
    ###################################################################################
    if plate==True and 'state' not in line.lower():
        cerberus=re.search('Dimensions ntax\=([0-9]+)\;',line) ## Extract useful information from the bits preceding the actual trees.
        if cerberus is not None:
            tipNum=int(cerberus.group(1))

        if 'Translate' in line:
            taxonlist=True ## taxon list to follow

        if taxonlist==True and ';' not in line and 'Translate' not in line: ## remember tip encodings
            cerberus=re.search('([0-9]+) ([\'\"A-Za-z0-9\?\|\-\_\.\/]+)',line)
            tips[cerberus.group(1)]=cerberus.group(2).strip("'")

    if 'tree STATE_' in line and plate==True: ## starting actual analysis
        plate=False
        assert (tipNum == len(tips)),'Expected number of tips: %s\nNumber of tips found: %s'%(tipNum,len(tips)) ## check that correct numbers of tips have been parsed
    ################################################################################### start analysing trees
    cerberus=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line) ## search for crud at the beginning of the line that's not a tree string

    if cerberus is not None: ## tree identified
        ################################################################# at state 0 - create the header for the output file and read the tree (in case the output log file requires information encoded in the tree)
        if treecount==0: ## At tree state 0 insert header into output file
            ll=bt.tree() ## empty tree object
            start=len(cerberus.group()) ## index of where tree string starts in the line
            treestring=str(line[start:]) ## grab tree string
            bt.make_tree(treestring,ll) ## read tree string
            if lower==0 and upper==np.inf: ## only add a header if not doing a chunk
                outfile.write('state') ## begin the output log file
                ########################################### add header to output log file
                if 'treeLength' in analyses:
                    outfile.write('\ttreeLength')
                ###########################################
                if 'RC' in analyses:
                    outfile.write('\tN\tS\tuN\tuS\tdNdS')
                ###########################################
                if 'tmrcas' in analyses:
                    tmrcas={'A':[],'B':[],'C':[]} ## dict of clade names
                    ll.renameTips(tips)
                    for k in ll.Objects: ## iterate over branches
                        if isinstance(k,bt.leaf): ## only interested in tips
                            if 'A' in k.name: ## if name of tip satisfies condition
                                tmrcas['A'].append(k.numName) ## add the tip's numName to tmrca list - tips will be used to ID the common ancestor
                            elif k.name in Conakry_tips:
                                tmrcas['B'].append(k.numName)
                            tmrcas['C'].append(k.numName)
                            
                    outfile.write('\t%s'%('\t'.join(sorted(tmrcas.keys()))))
                ###########################################
                if 'transitions' in analyses:
                    outfile.write('state\ttotalChangeCount\tcompleteHistory_1')
                ###########################################
                if 'subtrees' in analyses:
                    import copy
                ###########################################
                if 'subrate' in analyses:
                    subgroups={x:[] for x in substr} ## dict of leafs for subrate groups
                    subgroup_eye={'SLE_pC_eye_201403522_2014-12-14':'1'}
                    
                    for key, value in subgroup_eye.iteritems():
                        print key, value

                    #***************added below from Jason for info on monophyly*******************************
                    sg_mono={x:[] for x in substr} ## dict that will hold info about monophyly for the different groups in each of the trees
                    sg_leaves={x:{} for x in substr} ## dict that will hold info about the leaves of the MRCA node for each subgroup
                    #***************added above from Jason for info on monophyly*******************************
                    
                    ll.renameTips(tips)
                    for k in ll.Objects: ## iterate over branches
                        if isinstance(k,bt.leaf): ## only interested in tips
                            for s in subgroups.keys():
                                regexMatch=re.search(s, k.name)
                                if regexMatch:
                                    print regexMatch.group()
                                    print "found string"
                                    print k.name
                                    subgroups[s].append(k.numName)
                                    print subgroups[s]
                                    
                    #Find the branch for the eye rate.
                    for k in ll.Objects: ## iterate over branches
                        if isinstance(k,bt.leaf): ## only interested in tips
                            for s in subgroup_eye.keys():
                                regexMatch=re.search(s, k.name)
                                if regexMatch:
                                    print regexMatch.group()
                                    print "found eye string"
                                    print k.name
                                    subgroup_eye[s] = k.numName
                                    print subgroup_eye[s]
                                #if s in k.name: ## if name of tip satisfies condition
                                #     ## add the tip's numName to subrate group list - tips will be used to ID the common ancestor
                    outfile.write('\t%s\tEye_branch_rate\tBlood_Combo\tSemen_Combo\tAcute_Rate\tFull'%('\t'.join(sorted(subgroups.keys()))))

                ###########################################
                ## your custom header making code goes here
                ## if 'custom' in analyses:
                ##     trait='yourTrait'
                ##     trait_vals=[]
                ##     for k in ll.Objects:
                ##         if k.traits.has_key(trait):
                ##             trait_vals.append(k.traits[trait])
                ##     available_trait_values=sorted(bt.unique(trait_vals))
                ##     for tr in available_trait_values:
                ##         outfile.write('\t%s.time'%(tr))
                ###########################################
                treecount+1
                outfile.write('\n') ## newline for first tree
        #################################################################
        if int(cerberus.group(1)) >= burnin and lower <= int(cerberus.group(1)) < upper: ## After burnin start processing
            ll=bt.tree() ## ll is the tree object
            start=len(cerberus.group()) ## find start of tree string in line
            treestring=str(line[start:]) ## get tree string
            bt.make_tree(treestring,ll) ## pass it to make_tree function
            ll.traverse_tree() ## Traverse the tree - sets the height of each object in the tree
            #### renaming tips
            if len(tips)>0:
                ll.renameTips(tips) ## Rename tips so their name refers to sequence name
            else:
                for k in ll.Objects:
                    if isinstance(k,leaf):
                        k.name=k.numName ## otherwise every tip gets a name that's the same as tree string names
            #### calibration
            dateCerberus=re.compile('\|([0-9\-]+)$') ## search pattern + brackets on actual calendar date
            #print dateCerberus.groups
            if calibration==True: ## Calibrate tree so everything has a known position in actual time
                tipDatesRaw=[dateCerberus.search(x).group(1) for x in tips.values()]
                tipDates=map(bt.decimalDate,tipDatesRaw)
                maxDate=max(tipDates) ## identify most recent tip
                ll.setAbsoluteTime(maxDate)
            outfile.write('%s'%cerberus.group(1)) ## write MCMC state number to output log file
            ################################################################################
            if 'treeLength' in analyses:
                treeL=sum([k.length for k in ll.Objects]) ## do analysis
                outfile.write('\t%s'%(treeL)) ## output to file
            ###################################################
            if 'RC' in analyses: ## 'RC' was queued as an analysis
                Ns=[] ## empty list
                Ss=[]
                uNs=[]
                uSs=[]
                for k in ll.Objects: ## iterate over branch objects in the tree
                    if k.traits.has_key('N'): ## if branch has a trait labelled "N"...
                        Ns.append(k.traits['N']) ## add it to empty list
                        Ss.append(k.traits['S']) ## likewise for every other trait
                        uNs.append(k.traits['b_u_N'])
                        uSs.append(k.traits['b_u_S'])
                tNs=sum(Ns) ## sum of numbers in list
                tSs=sum(Ss)
                tuNs=sum(uNs)
                tuSs=sum(uSs)
                dNdS=(tNs/tSs)/(tuNs/tuSs) ## calculate dNdS
                outfile.write('\t%s\t%s\t%s\t%s\t%s'%(tNs,tSs,tuNs,tuSs,dNdS)) ## output to file, separated by tabs
            ###################################################
            if 'tmrcas' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                nodes={x:None for x in tmrcas.keys()} ## each TMRCA will correspond to a single object
                score={x:len(ll.Objects)+1 for x in tmrcas.keys()} ## this will be used to score the common ancestor candidate
                for required in tmrcas.keys(): ## iterate over TMRCAs
                    searchNodes=sorted([nd for nd in ll.Objects if nd.branchType=='node' and len(nd.leaves)>=len(tmrcas[required])],key=lambda n:len(n.leaves)) ## common ancestor candidates must have at least as many descendants as the list of search tips
                    for k in searchNodes: ## iterate over candidates
                        common,queryLeft,targetLeft=overlap(k.leaves,tmrcas[required]) ## find how many query tips exist as descendants of candidate nodes
                        if len(targetLeft)==0 and len(queryLeft)<=score[required]: ## all of query tips must be descended from common ancestor, every extra descendant of common ancestor not in query contributes to a score
                            nodes[required]=k ## if score improved - assign new common ancestor
                            score[required]=len(queryLeft) ## score is extra descendants not in the list of known tips
                        
                outTMRCA=['%.6f'%(nodes[n].absoluteTime) for n in sorted(nodes.keys())] ## fetch absoluteTime of each identified common ancestor
                outfile.write('\t%s'%('\t'.join(outTMRCA)))
            ###################################################
            if 'Sharp' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                assert len(analyses)==1,'More that one analysis queued in addition to Sharp, which is inadvisable'
                outSharp=[]
                for k in ll.Objects:
                    if k.traits.has_key('N'):
                        N=k.traits['N']
                        S=k.traits['S']
                        halfBranch=k.length*0.5
                        if isinstance(k,node):
                            all_leaves=[tips[lf] for lf in k.leaves]
                            t=min(map(bt.decimalDate,[dateCerberus.search(x).group(1) for x in all_leaves]))-k.absoluteTime+halfBranch
                        else:
                            t=halfBranch
                    
                        outSharp.append('(%d,%d,%.4f)'%(N,S,t))
                outfile.write('\t%s'%('\t'.join(outSharp)))
            ###################################################
            if 'transitions' in analyses:
                assert calibration==True,'This analysis type requires time-calibrated trees'
                assert len(analyses)==1,'More that one analysis queued in addition to transitions, which is inadvisable'
                outTransitions=[]
                for k in ll.Objects:
                    if k.traits.has_key('location.states') and k.parent.traits.has_key('location.states'):
                        cur_value=k.traits['location.states']
                        par_value=k.parent.traits['location.states']
                        if cur_value!=par_value:
                            outTransitions.append('{1,%s,%s,%s}'%(ll.treeHeight-k.height-0.5*k.length,par_value,cur_value))
                outfile.write('\t%d\t%s'%(len(outTransitions),'\t'.join(outTransitions)))
            ###################################################
            if 'subtrees' in analyses:
                traitName='location.states'
                assert [k.traits.has_key(traitName) for k in ll.Objects].count(True)>0,'No branches have the trait "%s"'%(traitName)
                for k in ll.Objects:
                    if k.traits.has_key(traitName) and k.parent.traits.has_key(traitName) and k.parent.index!='Root' and k.traits[traitName]!=k.parent.traits[traitName] and k.traits[traitName]=='human':
                        proceed=False ## assume we still can't proceed forward
                        kloc=k.traits[traitName]
                        if isinstance(k,bt.leaf): ## if dealing with a leaf - proceed
                            N_children=1
                            proceed=True
                        else:
                            N_children=len(k.leaves)
                            if [ch.traits[traitName] for ch in k.children].count(kloc)>0:
                                proceed=True
                        #print k.index,k.parent.index,k.traits,k.parent.traits,k.traits[traitName]
                    
                        if proceed==True: ## if at least one valid tip and no hanging nodes
                            subtree=copy.deepcopy(ll.traverseWithinTrait(k,traitName))
                            subtree_leaves=[x.name for x in subtree if isinstance(x,bt.leaf)]
                        
                            if len(subtree_leaves)>0:
                                mostRecentTip=max([bt.decimalDate(x.strip("'").split('|')[-1]) for x in subtree_leaves])
                                while sum([len(nd.children)-sum([1 if ch in subtree else 0 for ch in nd.children]) for nd in subtree if isinstance(nd,bt.node) and nd.index!='Root'])>0: ## keep removing nodes as long as there are nodes with children that are not entirely within subtree
                                    for nd in sorted([q for q in subtree if isinstance(q,bt.node)],key=lambda x:(sum([1 if ch in subtree else 0 for ch in x.children]),x.height)): ## iterate over nodes in subtree, starting with ones that have fewest valid children and are more recent
            
                                        child_status=[1 if ch in subtree else 0 for ch in nd.children] ## check how many children of current node are under the right trait value
            
                                        if sum(child_status)<2 and nd.index!='Root': ## if less than 2 children in subtree (i.e. not all children are under the same trait state)
                                            #print 'removing: %d, children in: %s'%(nd.index,[location_to_country[ch.traits[traitName]] for ch in nd.children])
                                            grand_parent=nd.parent ## fetch grandparent of node to be removed
                                            grand_parent.children.remove(nd) ## remove node from its parent's children

                                            if sum(child_status)==0: ## node has no valid children - current grandparent will be removed on next iteration, since it will only have one child
                                                pass
                                            else: ## at least one child is still valid - reconnect the one valid child to grandparent
                                                child=nd.children[child_status.index(1)] ## identify the valid child
                                                child.parent=grand_parent ## child's parent is now its grandparent
                                                grand_parent.children.append(child) ## child is now child of grandparent
                                                child.length+=nd.length ## child's length now includes it's former parent's length
                                            subtree.remove(nd) ## remove node from subtree
                                outfile.write('\t{%s,%s,%s,%s,%d}'%(k.absoluteTime,mostRecentTip,k.parent.traits[traitName],k.traits[traitName],len(subtree_leaves)))
                                sys.stderr.write('\t{%s,%s,%s,%s,%d}'%(k.absoluteTime,mostRecentTip,k.parent.traits[traitName],k.traits[traitName],len(subtree_leaves)))
                                ##########
                                ## Comment out to output stats rather than trees
                                ##########
#                                 if len(subtree)>0: ## only proceed if there's at least one tip in the subtree
#                                     local_tree=bt.tree() ## create a new tree object where the subtree will be
#                                     local_tree.Objects=subtree ## assign branches to new tree object
#                                     local_tree.root.children.append(subtree[0]) ## connect tree object's root with subtree
#                                     subtree[0].parent=local_tree.root ## subtree's root's parent is tree object's root
#                                     #local_tree.root.absoluteTime=subtree[0].absoluteTime-subtree[0].length ## root's absolute time is subtree's root time
#                                     local_tree.sortBranches() ## sort branches, draw small tree
#                                     subtreeString=local_tree.toString()
#                                     outfile.write('\t%s'%(subtreeString))
            ###################################################
            if 'subrate' in analyses:
                #assert calibration==True,'This analysis type requires time-calibrated trees'
                nodes=collections.OrderedDict()
                score=collections.OrderedDict()
                nodes={x:None for x in subgroups.keys()} ## each group will correspond to a single object
                score={x:len(ll.Objects)+1 for x in subgroups.keys()} ## this will be used to score the common ancestor candidate
                #This portion finds tmrca node for each subgroup
                for required in subgroups.keys(): ## iterate over groups
                    searchNodes=sorted([nd for nd in ll.Objects if nd.branchType=='node' and len(nd.leaves)>=len(subgroups[required])],key=lambda n:len(n.leaves)) ## common ancestor candidates must have at least as many descendants as the list of search tips
                    for k in searchNodes: ## iterate over candidates
                        common,queryLeft,targetLeft=overlap(k.leaves,subgroups[required]) ## find how many query tips exist as descendants of candidate nodes
                        if len(targetLeft)==0 and len(queryLeft)<=score[required]: ## all of query tips must be descended from common ancestor, every extra descendant of common ancestor not in query contributes to a score
                            nodes[required]=k ## if score improved - assign new common ancestor
                            score[required]=len(queryLeft) ## score is extra descendants not in the list of known tips
                            
                #Print out info about # taxa in MRCA clade for each group
                #for g,s in score.iteritems(): ## iterate over groups, nodes
                #    if s == 0: sg_mono[g].append(float(1))
                #    else: 
                #        sg_mono[g].append(float(0))
                #    clade_leaves = tuple(sorted([tips[l] for l in nodes[g].leaves]))
                #    if clade_leaves not in sg_leaves[g]: sg_leaves[g][clade_leaves]=1
                #    else: sg_leaves[g][clade_leaves] += 1
                
                #Old length and substitutions counting function is below: 
                #NOTE: it calculates rates from all branches of tMRCA, even is not monophyletic.
                #lensubs={x:None for x in subgroups.keys()} ## total branch lengths, substitutions for subgroups
                #for g,n in nodes.iteritems(): ## iterate over groups, nodes
                #    lensubs[g] = subclade_branch_info(n)
#                    print treecount, lensubs[g]

                #Combine the nodes and subgroups dictionaries by keys.
                #nodes_and_leaves dictionary will be passed to lineage_trace_info.
                #nodes_and_leaves = {'LBR_pA_plasma': [node_TMRCA],['103', '104'],
                #                    'SLE_pC_blood': [node_TMRCA],['372', '369', '368', '367', '371', '370'],
                #                    'GIN_pE_blood|GIN_pE_plasma': [node_TMRCA],['55', '57', '56'] }
                
                
                nodes_and_leaves=collections.OrderedDict()
                for key in (nodes.viewkeys() | subgroups.keys()):
                    if key in nodes: nodes_and_leaves.setdefault(key, []).append(nodes[key])
                    if key in subgroups: nodes_and_leaves.setdefault(key, []).append(subgroups[key])
                
                
                #Initialize the dictionary with subgroup.key values - i.e. user supplied search strings
                #and empty values.
                branches_want=collections.OrderedDict()
                branches_want={x:None for x in subgroups.keys()} ## total branch lengths, substitutions for subgroups - only leaves we want.
                #Print the contents of the nodes_and_leaves dictionary:
                for g,n in nodes_and_leaves.iteritems():
                    #print g     #user-supplied keys.
                    #print n[0]  #tMRCA_node
                    #print n[1]  #leaves_we_want_to_keep
                    
                    #branches_want[g] = lineage_trace_info(tMRCA_node, leaves_we_want)
                    ##***************Need to submit the blood and semen sequence list in interspersed order!!!**********
                    #Submit as: pA_blood pA_semen pC_blood pC_semen pE_blood pE_semen
                    #
                    #
                    #branches_want[g] = {'LBR_pA_plasma': [list of branches from tMRCA_pA_blood],
                    #                    'pA_semen': [list of branches from tMRCA_pA_semen],
                    #                    'SLE_pC_blood': [list of branches from tMRCA_pC_blood],
                    #                    'pC_semen': [list of branches from tMRCA_pC_semen]
                    #                    'GIN_pE_blood|GIN_pE_plasma': [list of branches from tMRCA_pE_blood],
                    #                    'pE_semen': [list of branches from tMRCA_pE_blood] }
                    #NOTE: Data is stored in orderedDict(branches_want) as FILO.  pA_blood pA_semen
                    #This becomes: 
                    #pA_semen 0
                    #pA_blood 1
                    #*****************************************************************************************************
                    
                    lineage_trace = []
                    branches_want[g] = lineage_trace_info(n[0], n[1])
                
                #The for loop below trims the intermingled branches between sets.
                #We place the shared branches into the blood set.  
                #IE - This is an assumption that shared blood and semen mutations originally arose
                #from the blood during acute infection.
                
                #For loop for troubleshooting:
                #for g,n in branches_want.iteritems():
                    #print "Dictionary key order is: %s" %(g)
                    #print "Dictionary value order is: %s" %(n)
                branches_want_ordered=collections.OrderedDict(sorted(branches_want.items(), reverse=True))
                #print branches_want
                #print branches_want_ordered
                
                for g, i in enumerate(branches_want_ordered) :
                    while (( g ) < len(branches_want_ordered) and (g % 2 == 0) ):
                        #Compare blood and semen lineage_trace sets.
                        #Remove shared branches from semen lineage_trace.
                        #n[g]=blood
                        #n[g+1]=semen
                        
                        #print g
                        #print i
                        #print branches_want.values()[g]
                        #print "Branches want[%d] %s" %(g, branches_want_ordered.values()[g])
                        #print "Branches want[%d] %s" %((g+1), branches_want_ordered.values()[g+1])
                        
                        #This is where the magic happens:
                        #Branches shared between semen and blood sets are pruned from semen sets and kept in blood set.
                        branches_want_ordered.values()[g].difference_update(branches_want_ordered.values()[g+1])
                    
                        #print "Updated branch[%d] %s" %(g, branches_want_ordered.values()[g])
                        g+=1
                        #print g
                    
                lensubs_want=collections.OrderedDict()
                lensubs_want={x:None for x in subgroups.keys()} ## total branch lengths, substitutions for subgroups - only leaves we want.
                #Calculate evolutionary rates trimmed branches:.
                i = 0
                for r,s in branches_want_ordered.iteritems():
                    #Initialize variables.
                    lengths=0
                    numsubs=0
                    #print r
                    #print s

                    for lin in s:
                        #print lin
                        lengths+=lin.length #lengths
                        numsubs+=lin.length * lin.traits['rate']
                    #Add an iterator value to the ordered dictionary to slice out specific blood and semen entries.
                    lensubs_want[r]=lengths,numsubs, i
                    i+=1
                    
                    #print "G is: %s" %(r)
                    #print lensubs_want[r][1]
                    #print lensubs_want[r][0]
                    #print lensubs_want[r][2]
                    #print "Evolutionary rate is: %s" %(lensubs_want[r][1] / lensubs_want[r][0])
                    

                #outSubrate=['%.6f'%(lensubs[x][1]/lensubs[x][0]) for x in sorted(nodes.keys())] ## Calculate average rates
                outSubrate=['%.6f'%(lensubs_want[x][1]/lensubs_want[x][0]) for x in sorted(nodes.keys())] ## Calculate average rates
                
                #Copy contents of lensubs_want to blood and semen-specific dictionaries
                #to calculate blood and semen-specific rates.
                lensubs_want_blood={}
                lensubs_want_semen={}
                lensubs_want_blood = {iterator:iterator_value for iterator, iterator_value in lensubs_want.iteritems() if ((iterator_value[2] % 2) == 1)} 
                #print lensubs_want_blood
                lensubs_want_semen = {iterator:iterator_value for iterator, iterator_value in lensubs_want.iteritems() if ((iterator_value[2] % 2) == 0)}
                #print lensubs_want_semen
                    
                #Modified below to calculate blood- and semen-specific combo rates.
                combo_blood_SubgroupsRate = sum([x[1] for x in lensubs_want_blood.values()])/sum([x[0] for x in lensubs_want_blood.values()])
                #print "blood_SubgroupsRate: %.6f " %(combo_blood_SubgroupsRate)
                combo_semen_SubgroupsRate = sum([x[1] for x in lensubs_want_semen.values()])/sum([x[0] for x in lensubs_want_semen.values()])
                #print "semen_SubgroupsRate: %.6f " %(combo_semen_SubgroupsRate)
                
                #Calculate the rate for the branch leading to the eye sample.
                eye_branch_rate_values={x:None for x in subgroup_eye.keys()} ## total branch lengths, substitutions for subgroups - only leaves we want.
                for eye_key, eye_value in subgroup_eye.iteritems():
                    ll.renameTips(tips)
                    for k in ll.Objects: #Iterate over branches
                        if isinstance(k, bt.leaf):
                            if eye_value == k.numName:
                                #print "Found eye branch: %s" %(eye_value)
                                eye_branch_rate_values[eye_key] = k.length, (k.length * k.traits['rate'])
                    
                    #print "Eye substitution rates:"
                    #print eye_branch_rate_values[eye_key][0]
                    #print eye_branch_rate_values[eye_key][1]
                    eye_branch_rate = eye_branch_rate_values[eye_key][1] / eye_branch_rate_values[eye_key][0]
                    #print "Eye rate is: %.6f" %(eye_branch_rate)
                    
                    #Add the eye values to the lensubs_want dictionary.
                    #combo_SubgroupsRate holds all of the persistence rates!
                    #print "Adding value to lensubs_want_semen"
                    lensubs_want_semen[eye_key] = eye_branch_rate_values[eye_key][0], eye_branch_rate_values[eye_key][1]
                    #print lensubs_want_semen
                
                combo_Persistence_SubgroupsRate = sum([x[1] for x in lensubs_want_semen.values()])/sum([x[0] for x in lensubs_want_semen.values()])


                #Get info for full tree
                fullLen, fullSubs = subclade_branch_info(ll.root)
                #Calc Avg across "Rest" of tree (i.e., all branches not included in any subgroup)
                #!!!!! Assumes that there is no overlap between subgroups
                restRate = (fullSubs - sum([x[1] for x in lensubs_want_semen.values()]))/(fullLen - sum([x[0] for x in lensubs_want_semen.values()]))
                outfile.write('\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f'%('\t'.join(outSubrate), eye_branch_rate, combo_blood_SubgroupsRate, combo_semen_SubgroupsRate, restRate, fullSubs/fullLen))
            ###################################################
            ## your analysis and output code goes here, e.g.
            ## if 'custom' in analyses:
            ##     out={x:0.0 for x in available_trait_values}
            ##     for k in ll.Objects:
            ##         if k.traits.has_key(trait):
            ##             out[k.traits[trait]]+=k.length
            ##     for tr in available_trait_values:
            ##         outfile.write('\t%s'%(out[tr]))
            ###################################################
            outfile.write('\n') ## newline for post-burnin tree
        treecount+=1 ## increment tree counter
        ################################################################################            
        if treecount==threshold: ## tree passed progress bar threshold
            timeTakenSoFar=dt.datetime.now()-begin ## time elapsed
            timeElapsed=float(divmod(timeTakenSoFar.total_seconds(),60)[0]+(divmod(timeTakenSoFar.total_seconds(),60)[1])/float(60))
            timeRate=float(divmod(timeTakenSoFar.total_seconds(),60)[0]*60+divmod(timeTakenSoFar.total_seconds(),60)[1])/float(treecount+1) ## rate at which trees have been processed
            processingRate.append(timeRate) ## remember rate
            ETA=(sum(processingRate)/float(len(processingRate))*(Ntrees-treecount))/float(60)/float(60) ## estimate how long it'll take, given mean processing rate
    
            excessiveTrees=treecount
            if treecount>=10000:
                excessiveTrees=10000
            if timeElapsed>60.0: ## took over 60 minutes
                reportElapsed=timeElapsed/60.0 ## switch to hours
                reportUnit='h' ## unit is hours
            else:
                reportElapsed=timeElapsed ## keep minutes
                reportUnit='m'
        
            sys.stderr.write('\r') ## output progress bar
            sys.stderr.write("[%-30s] %4d%%  trees: %5d  elapsed: %5.2f%1s  ETA: %5.2fh (%6.1e s/tree)" % ('='*(excessiveTrees/progress_update),treecount/float(Ntrees)*100.0,treecount,reportElapsed,reportUnit,ETA,processingRate[-1]))
            sys.stderr.flush()

            threshold+=progress_update ## increment to next threshold
            
            #*****************************Added below from Jason's script*********************************************
            #if 'subrate' in analyses:
            #    print ("\n%s" %"\t".join(sg_mono.keys()))
            #    print "\t".join(["%.2f" % (np.mean(x)*100) for x in sg_mono.values()])
            #    for g, ldict in sg_leaves.iteritems():
            #        print "\n%s:" % g
            #        for k, v in ldict.iteritems():
            #            print "%s\t%.2f" % (k, (float(v)/len(sg_mono[g]))*100)
            #*****************************Added above from Jason's script*********************************************
            
        ################################################################################
        if 'End;' in line:
            pass
outfile.close()
sys.stderr.write('\nDone!\n') ## done!

#############################################
