#!/usr/bin/env python


# Code based and adapted from https://github.com/ljdursi/poapy

from __future__ import print_function
try:
    from builtins import zip
    from builtins import str
    from builtins import object
    from builtins import range
except ImportError:
    pass
import numpy
import textwrap
import collections


class Node(object):
    def __init__(self, nodeID=-1, base='N'):
        self.ID = nodeID
        self.base = base
        self.inEdges = {}
        self.outEdges = {}
        self.alignedTo = []

    def __str__(self):
        return "(%d:%s)" % (self.ID, self.base)

    def _add_edge(self, edgeset, neighbourID, label, from_neighbour):
        if neighbourID is None:
            return
        # already present? just update labels
        # otherwise create appropriately-ordered edge and proceed
        if neighbourID in edgeset:
            edgeset[neighbourID].addLabel(label)
        else:
            if from_neighbour:
                edge = Edge(outNodeID=neighbourID, inNodeID=self.ID, label=label)
            else:
                edge = Edge(outNodeID=self.ID, inNodeID=neighbourID, label=label)
            edgeset[neighbourID] = edge

    def addInEdge(self, neighbourID, label):
        self._add_edge(self.inEdges, neighbourID, label, from_neighbour=True)

    def addOutEdge(self, neighbourID, label):
        self._add_edge(self.outEdges, neighbourID, label, from_neighbour=False)

    def nextNode(self, label):
        """Returns the first (presumably only) outward neighbour
           having the given edge label"""
        nextID = None
        for e in self.outEdges:
            if label in self.outEdges[e].labels:
                nextID = e
        return nextID

    @property
    def inDegree(self):
        return len(self.inEdges)

    @property
    def outDegree(self):
        return len(self.outEdges)

    @property
    def labels(self):
        """Returns all the labels associated with an in-edge or an out edge."""
        labelset = set([])
        for e in list(self.inEdges.values()):
            labelset = labelset.union(e.labels)
        for e in list(self.outEdges.values()):
            labelset = labelset.union(e.labels)
        return list(labelset)


class Edge(object):
    def __init__(self, inNodeID=-1, outNodeID=-1, label=None):
        self.inNodeID  = inNodeID
        self.outNodeID = outNodeID
        if label is None:
            self.labels = []
        elif type(label) == list:
            self.labels = label
        else:
            self.labels = [label]

    def addLabel(self, newlabel):
        if newlabel not in self.labels:
            self.labels.append(newlabel)
        return

    def __str__(self):
        nodestr = "(%d) -> (%d) " % (self.outNodeID, self.inNodeID)
        if self.labels is None:
            return nodestr
        else:
            return nodestr + self.labels.__str__()


class POAGraph(object):
    def addUnmatchedSeq(self, seq, label=None, updateSequences=True):
        """Add a completely independant (sub)string to the graph,
           and return node index to initial and final node"""
        if seq is None:
            return

        firstID, lastID = None, None
        neededSort = self.needsSort

        for base in seq:
            nodeID = self.addNode(base)
            if firstID is None:
                firstID = nodeID
            if lastID is not None:
                self.addEdge(lastID, nodeID, label)
            lastID = nodeID

        self.__needsort = neededSort  # no new order problems introduced
        if updateSequences:
            self.__seqs.append(seq)
            self.__labels.append(label)
            self.__starts.append(firstID)
        return firstID, lastID

    def __init__(self, seq=None, label=None):
        self._nextnodeID = 0
        self._nnodes = 0
        self._nedges = 0
        self.nodedict = {}
        self.nodeidlist = []   # allows a (partial) order to be imposed on the nodes
        self.__needsort = False
        self.__labels = []
        self.__seqs = []
        self.__starts = []

        if seq is not None:
            self.addUnmatchedSeq(seq, label)

    def nodeIdxToBase(self, idx):
        return self.nodedict[self.nodeidlist[idx]].base

    def addNode(self, base):
        nid = self._nextnodeID
        newnode = Node(nid, base)
        self.nodedict[nid] = newnode
        self.nodeidlist.append(nid)
        self._nnodes += 1
        self._nextnodeID += 1
        self._needsSort = True
        return nid

    def addEdge(self, start, end, label):
        if start is None or end is None:
            return

        if start not in self.nodedict:
            raise KeyError('addEdge: Start node not in graph: '+str(start))
        if end not in self.nodedict:
            raise KeyError('addEdge: End node not in graph: '+str(end))

        oldNodeEdges = self.nodedict[start].outDegree + self.nodedict[end].inDegree

        self.nodedict[start].addOutEdge(end, label)
        self.nodedict[end].addInEdge(start, label)

        newNodeEdges = self.nodedict[start].outDegree + self.nodedict[end].inDegree

        if newNodeEdges != oldNodeEdges:
            self._nedges += 1

        self._needsSort = True
        return

    @property
    def needsSort(self):
        return self.__needsort

    @property
    def nNodes(self):
        return self._nnodes

    @property
    def nEdges(self):
        return self._nedges

    def _simplified_graph_rep(self):
        ## TODO: The need for this suggests that the way the graph is currently represented
        ## isn't really right and needs some rethinking.

        node_to_pn = {}
        pn_to_nodes = {}

        # Find the mappings from nodes to pseudonodes
        cur_pnid = 0
        for _, node in self.nodedict.items():
            if node.ID not in node_to_pn:
                node_ids = [node.ID] + node.alignedTo
                pn_to_nodes[cur_pnid] = node_ids
                for nid in node_ids:
                    node_to_pn[nid] = cur_pnid
                cur_pnid += 1

        # create the pseudonodes
        Pseudonode = collections.namedtuple("Pseudonode", ["pnode_id", "predecessors", "successors", "node_ids"])
        pseudonodes = []

        for pnid in range(cur_pnid):
            nids, preds, succs = pn_to_nodes[pnid], [], []
            for nid in nids:
                node = self.nodedict[nid]
                preds += [node_to_pn[inEdge.outNodeID] for _, inEdge in node.inEdges.items()]
                succs += [node_to_pn[outEdge.inNodeID] for _, outEdge in node.outEdges.items()]

            pn = Pseudonode(pnode_id=pnid, predecessors=preds, successors=succs, node_ids=nids)
            pseudonodes.append(pn)

        return pseudonodes

    def toposort(self):
        """Sorts node list so that all incoming edges come from nodes earlier in the list."""
        sortedlist = []
        completed = set([])

        ##
        ## The topological sort of this graph is complicated by the alignedTo edges;
        ## we want to nodes connected by such edges to remain near each other in the
        ## topological sort.
        ##
        ## Here we'll create a simple version of the graph that merges nodes that
        ## are alignedTo each other, performs the sort, and then decomposes the
        ## 'pseudonodes'.
        ##
        ## The need for this suggests that the way the graph is currently represented
        ## isn't quite right and needs some rethinking.
        ##

        pseudonodes = self._simplified_graph_rep()

        def dfs(start, complete, sortedlist):
            stack, started = [start], set()
            while stack:
                pnodeID = stack.pop()

                if pnodeID in complete:
                    continue

                if pnodeID in started:
                    complete.add(pnodeID)
                    for nid in pseudonodes[pnodeID].node_ids:
                        sortedlist.insert(0, nid)
                    started.remove(pnodeID)
                    continue

                successors = pseudonodes[pnodeID].successors
                started.add(pnodeID)
                stack.append(pnodeID)
                stack.extend(successors)

        while len(sortedlist) < self.nNodes:
            found = None
            for pnid in range(len(pseudonodes)):
                if pnid not in completed and len(pseudonodes[pnid].predecessors) == 0:
                    found = pnid
                    break
            assert found is not None
            dfs(found, completed, sortedlist)

        assert len(sortedlist) == self.nNodes
        self.nodeidlist = sortedlist
        self._needsSort = False
        return

    def testsort(self):
        """ Test the nodeidlist to make sure it is topologically sorted:
            eg, all predecessors of a node preceed the node in the list"""
        if self.nodeidlist is None:
            return
        seen_nodes = set()
        for nodeidx in self.nodeidlist:
            node = self.nodedict[nodeidx]
            for in_neighbour in node.inEdges:
                assert in_neighbour in seen_nodes
            seen_nodes.add(nodeidx)
        return

    def nodeiterator(self):
        if self.needsSort:
            self.toposort()

        def nodegenerator():
            for nodeidx in self.nodeidlist:
                yield self.nodedict[nodeidx]

        return nodegenerator

    def __str__(self):
        selfstr = ""
        ni = self.nodeiterator()
        for node in ni():
            selfstr += node.__str__() + "\n"
            for outIdx in node.outEdges:
                selfstr += "        " + node.outEdges[outIdx].__str__() + "\n"
        return selfstr

    def incorporateSeqAlignment(self, alignment, seq, label=None):
        """Incorporate a SeqGraphAlignment into the graph."""
        newseq     = alignment.sequence
        stringidxs = alignment.stringidxs
        nodeidxs   = alignment.nodeidxs

        firstID = None
        headID = None
        tailID = None

        # head, tail of sequence may be unaligned; just add those into the
        # graph directly
        validstringidxs = [si for si in stringidxs if si is not None]
        startSeqIdx, endSeqIdx = validstringidxs[0], validstringidxs[-1]
        if startSeqIdx > 0:
            firstID, headID = self.addUnmatchedSeq(newseq[0:startSeqIdx], label, updateSequences=False)
        if endSeqIdx < len(newseq):
            tailID, __ = self.addUnmatchedSeq(newseq[endSeqIdx+1:], label, updateSequences=False)

        # now we march along the aligned part. For each base, we find or create
        # a node in the graph:
        #   - if unmatched, the corresponding node is a new node
        #   - if matched:
        #       - if matched to a node with the same base, the node is that node
        #       - if matched to a node with a different base whch is in turn
        #         aligned to a node with the same base, that aligned node is
        #         the node
        #       - otherwise, we create a new node.
        # In all cases, we create edges (or add labels) threading through the
        # nodes.
        for sindex, matchID in zip(stringidxs, nodeidxs):
            if sindex is None:
                continue
            base = newseq[sindex]
            if matchID is None:
                nodeID = self.addNode(base)
            elif self.nodedict[matchID].base == base:
                nodeID = matchID
            else:
                otherAligns = self.nodedict[matchID].alignedTo
                foundNode = None
                for otherNodeID in otherAligns:
                    if self.nodedict[otherNodeID].base == base:
                        foundNode = otherNodeID
                if foundNode is None:
                    nodeID = self.addNode(base)
                    self.nodedict[nodeID].alignedTo = [matchID] + otherAligns
                    for otherNodeID in [matchID] + otherAligns:
                        self.nodedict[otherNodeID].alignedTo.append(nodeID)
                else:
                    nodeID = foundNode

            self.addEdge(headID, nodeID, label)
            headID = nodeID
            if firstID is None:
                firstID = headID

        # finished the unaligned portion: now add an edge from the current headID to the tailID.
        self.addEdge(headID, tailID, label)

        # resort
        self.toposort()

        self.__seqs.append(seq)
        self.__labels.append(label)
        self.__starts.append(firstID)
        return

    def consensus(self, excludeLabels=None):
        if excludeLabels is None:
            excludeLabels = []

        if self.needsSort:
            self.toposort()

        nodesInReverse = self.nodeidlist[::-1]
        maxnodeID = max(nodesInReverse)+1
        nextInPath = [-1]*maxnodeID
        scores = numpy.zeros((maxnodeID))

        for nodeID in nodesInReverse:
            bestWeightScoreEdge = (-1, -1, None)

            for neighbourID in self.nodedict[nodeID].outEdges:
                e = self.nodedict[nodeID].outEdges[neighbourID]
                weight = len([l for l in e.labels if l not in excludeLabels])
                weightScoreEdge = (weight, scores[neighbourID], neighbourID)

                if weightScoreEdge > bestWeightScoreEdge:
                    bestWeightScoreEdge = weightScoreEdge

            scores[nodeID] = sum(bestWeightScoreEdge[0:2])
            nextInPath[nodeID] = bestWeightScoreEdge[2]

        pos = numpy.argmax(scores)
        path   = []
        bases  = []
        labels = []
        while pos is not None and pos > -1:
            path.append(pos)
            bases.append(self.nodedict[pos].base)
            labels.append(self.nodedict[pos].labels)
            pos = nextInPath[pos]

        return path, bases, labels

    def allConsenses(self, maxfraction=0.5):
        allpaths = []
        allbases = []
        alllabels = []
        exclusions = []

        passno = 0
        lastlen = 1000
        maxpasses = 10

        while len(exclusions) < len(self.__labels) and lastlen >= 10 and passno < maxpasses:
            path, bases, labellists = self.consensus(exclusions)
            if len(path) > 0:
                allpaths.append(path)
                allbases.append(bases)
                alllabels.append(labellists)

                labelcounts = collections.defaultdict(int)
                for ll in labellists:
                    for l in ll:
                        labelcounts[l] += 1

                for label, seq in zip(self.__labels, self.__seqs):
                    if label in labelcounts and labelcounts[
                        label
                    ] >= maxfraction * len(seq):
                        exclusions.append(label)

            lastlen = len(path)
            passno += 1

        return list(zip(allpaths, allbases, alllabels))

    def generateAlignmentStrings(self):
        """ Return a list of strings corresponding to the alignments in the graph """

        # Step 1: assign node IDs to columns in the output
        #  column_index[node.ID] is the position in the toposorted node list
        #    of the node itself, or the earliest node it is aligned to.
        column_index = {}
        current_column = 0

        # go through nodes in toposort order
        ni = self.nodeiterator()
        for node in ni():
            other_columns = [column_index[other] for other in node.alignedTo if other in column_index]
            if other_columns:
                found_idx = min(other_columns)
            else:
                found_idx = current_column
                current_column += 1

            column_index[node.ID] = found_idx

        ncolumns = current_column

        # Step 2: given the column indexes, populate the strings
        #   corresponding to the sequences inserted in the graph
        seqnames = []
        alignstrings = []
        for label, start in zip(self.__labels, self.__starts):
            seqnames.append(label)
            curnode_id = start
            charlist = ['  -  ']*ncolumns
            while curnode_id is not None:
                node = self.nodedict[curnode_id]
                if False:
                    print("str")
                    sep = ""
                    charlist[column_index[curnode_id]] = node.base
                else:
                    sep=","
                    charlist[column_index[curnode_id]] = f"{node.base:+5d}".replace('+', ' ')
                curnode_id = node.nextNode(label)
            alignstrings.append(sep.join(charlist))

        # Step 3: Same as step 2, but with consensus sequences
        consenses = self.allConsenses()
        for i, consensus in enumerate(consenses):
            seqnames.append('Consensus'+str(i))
            charlist = ['  -  ']*ncolumns
            for path, base in zip(consensus[0], consensus[1]):
                sep=","
                charlist[column_index[path]] = f"{base:+5d}".replace('+', ' ')
                # charlist[column_index[path]] = base
            alignstrings.append(sep.join(charlist))

        return list(zip(seqnames, alignstrings))
    
    def generateAlignment(self):
        """ Return a list of strings corresponding to the alignments in the graph """

        # Step 1: assign node IDs to columns in the output
        #  column_index[node.ID] is the position in the toposorted node list
        #    of the node itself, or the earliest node it is aligned to.
        column_index = {}
        current_column = 0

        # go through nodes in toposort order
        ni = self.nodeiterator()
        for node in ni():
            other_columns = [column_index[other] for other in node.alignedTo if other in column_index]
            if other_columns:
                found_idx = min(other_columns)
            else:
                found_idx = current_column
                current_column += 1

            column_index[node.ID] = found_idx

        ncolumns = current_column

        set_list = [set() for _ in range(ncolumns)]

        # Step 2: given the column indexes, populate the strings
        #   corresponding to the sequences inserted in the graph
        for label, start in zip(self.__labels, self.__starts):
            curnode_id = start
            num_list = [None]*ncolumns
            while curnode_id is not None:
                node = self.nodedict[curnode_id]
                num_list[column_index[curnode_id]]= node.base
                curnode_id = node.nextNode(label)
            for i, base in enumerate(num_list):
                set_list[i].add(base)
            
        return set_list

    def jsOutput(self):
        """returns a list of strings containing a a description of the graph for viz.js, http://visjs.org"""

        # get the consensus sequence, which we'll use as the "spine" of the
        # graph
        path, __, __ = self.consensus()
        pathdict = {nodeID: i*150 for i, nodeID in enumerate(path)}
        lines = ['var nodes = [']

        ni = self.nodeiterator()
        count = 0
        for node in ni():
            line = '    {id:'+str(node.ID)+', label: "'+str(node.base)+'"'
            if node.ID in pathdict and count % 5 == 0:
                line += ', x: ' + str(pathdict[node.ID]) + ', y: 0 , fixed: { x:true, y:false }},'
            else:
                line += '},'
            lines.append(line)

        lines[-1] = lines[-1][:-1]
        lines.append('];')

        lines.append(' ')

        lines.append('var edges = [')
        ni = self.nodeiterator()
        for node in ni():
            nodeID = str(node.ID)
            for edge in node.outEdges:
                target = str(edge)
                weight = str(len(node.outEdges[edge].labels)+1)
                lines.append('    {from: '+nodeID+', to: '+target+', value: '+weight+'},')
            for alignededge in node.alignedTo:
                # These edges indicate alignment to different bases, and are
                # undirected; thus make sure we only plot them once:
                if node.ID > alignededge:
                    continue
                target = str(alignededge)
                lines.append('    {from: '+nodeID+', to: '+target+', value: 1, color: "red", dashes: true, style: "dash-line"},')
        lines[-1] = lines[-1][:-1]
        lines.append('];')
        return lines

    def htmlOutput(self, outfile):
        header = """
                  <!doctype html>
                  <html>
                  <head>
                    <title>POA Graph Alignment</title>

                    <script type="text/javascript" src="https://unpkg.com/vis-network@9.0.4/standalone/umd/vis-network.min.js"></script>
                  </head>

                  <body>

                  <div id="loadingProgress">0%</div>

                  <div id="mynetwork"></div>

                  <script type="text/javascript">
                    // create a network
                  """
        outfile.write(textwrap.dedent(header[1:]))
        lines = self.jsOutput()
        for line in lines:
            outfile.write(line+'\n')
        footer = """
                  var container = document.getElementById('mynetwork');
                  var data= {
                    nodes: nodes,
                    edges: edges,
                  };
                  var options = {
                    width: '100%',
                    height: '800px',
                    physics: {
                        stabilization: {
                            updateInterval: 10,
                        }
                    }
                  };
                  var network = new vis.Network(container, data, options);

                  network.on("stabilizationProgress", function (params) {
                    document.getElementById("loadingProgress").innerText = Math.round(params.iterations / params.total * 100) + "%";
                  });
                  network.once("stabilizationIterationsDone", function () {
                      document.getElementById("loadingProgress").innerText = "100%";
                      setTimeout(function () {
                        document.getElementById("loadingProgress").style.display = "none";
                      }, 500);
                  });

                </script>

                </body>
                </html>
                """
        outfile.write(textwrap.dedent(footer))

class SeqGraphAlignment(object):
    __matchscore = 1
    __mismatchscore = -1
    __gap = -2

    def __init__(self, sequence, graph, fastMethod=True, globalAlign=False,
                 matchscore=__matchscore, mismatchscore=__mismatchscore,
                 gapscore=__gap, *args, **kwargs):
        self._mismatchscore = mismatchscore
        self._matchscore = matchscore
        self._gap = gapscore
        self.sequence    = sequence
        self.graph       = graph
        self.stringidxs  = None
        self.nodeidxs    = None
        self.globalAlign = globalAlign
        if fastMethod:
            matches = self.alignStringToGraphFast(*args, **kwargs)
        else:
            matches = self.alignStringToGraphSimple(*args, **kwargs)
        self.stringidxs, self.nodeidxs = matches

    def alignmentStrings(self):
        return "".join(
            str(self.sequence[i]) if i is not None else "  -  " for i in self.stringidxs
        ), "".join(
            self.graph.nodedict[j].base if j is not None else "  -  "
            for j in self.nodeidxs
        )

    def matchscore(self, c1, c2):
        if c1 == c2:
            return self._matchscore
        else:
            return self._mismatchscore

    def matchscoreVec(self, c, v):
        return numpy.where(v == c, self._matchscore, self._mismatchscore)

    def alignStringToGraphSimple(self):
        """Align string to graph, following same approach as smith waterman
        example"""
        # if type(self.sequence) != str:
        #     raise TypeError("Invalid Type")

        nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i, node in enumerate(ni()):
            pbase = node.base

            for j, sbase in enumerate(self.sequence):
                # add all candidates to a list, pick the best
                candidates = [(scores[i+1, j] + self._gap, i+1, j, "INS")]
                for predIndex in self.prevIndices(node, nodeIDtoIndex):
                    candidates += [(scores[predIndex+1, j+1] + self._gap, predIndex+1, j+1, "DEL")]
                    candidates += [(scores[predIndex+1, j] + self.matchscore(sbase, pbase), predIndex+1, j, "MATCH")]

                scores[i+1, j+1], backGrphIdx[i+1, j+1], backStrIdx[i+1, j+1], movetype = max(candidates)

                if not self.globalAlign and scores[i+1, j+1] < 0:
                    scores[i+1, j+1] = 0.
                    backGrphIdx[i+1, j+1] = -1
                    backStrIdx[i+1, j+1] = -1

        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)

    def alignStringToGraphFast(self):
        """Align string to graph - using numpy to vectorize across the string
        at each iteration."""
        # if not type(self.sequence) == str:
        #     raise TypeError("Invalid Type")

        l2 = len(self.sequence)
        seqvec = numpy.array(list(self.sequence))

        nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx = self.initializeDynamicProgrammingData()
        inserted = numpy.zeros((l2), dtype=numpy.bool)

        # having the inner loop as a function improves performance
        # can use Cython, etc on this for significant further improvements
        # can't vectorize this since there's a loop-carried dependency
        #  along the string
        def insertions(i, l2, scores, inserted):
            inserted[:] = False
            for j in range(l2):
                insscore = scores[i+1, j] + self._gap
                if insscore >= scores[i+1, j+1]:
                    scores[i+1, j+1] = insscore
                    inserted[j] = True

        # Dynamic Programming
        ni = self.graph.nodeiterator()
        for i, node in enumerate(ni()):
            gbase = node.base
            predecessors = self.prevIndices(node, nodeIDtoIndex)

            # calculate all best deletions, matches in one go over all
            # predecessors.

            # First calculate for the first predecessor, over all string posns:
            deletescore = scores[predecessors[0]+1, 1:] + self._gap
            bestdelete = numpy.zeros((l2), dtype=int)+predecessors[0]+1

            matchpoints = self.matchscoreVec(gbase, seqvec)
            matchscore = scores[predecessors[0]+1, 0:-1] + matchpoints
            bestmatch = numpy.zeros((l2), dtype=int)+predecessors[0]+1

            # then, the remaining
            for predecessor in predecessors[1:]:
                newdeletescore = scores[predecessor+1, 1:] + self._gap
                bestdelete     = numpy.where(newdeletescore > deletescore, predecessor+1, bestdelete)
                deletescore    = numpy.maximum(newdeletescore, deletescore)

                gbase = self.graph.nodeIdxToBase(predecessor)
                matchpoints = self.matchscoreVec(gbase, seqvec)
                newmatchscore = scores[predecessor+1, 0:-1] + matchpoints
                bestmatch     = numpy.where(newmatchscore > matchscore, predecessor+1, bestmatch)
                matchscore    = numpy.maximum(newmatchscore, matchscore)

            # choose best options available of match, delete
            deleted       = deletescore >= matchscore
            backGrphIdx[i+1, 1:] = numpy.where(deleted, bestdelete, bestmatch)
            backStrIdx [i+1, 1:] = numpy.where(deleted, numpy.arange(1, l2+1), numpy.arange(0, l2))
            scores[i+1, 1:] = numpy.where(deleted, deletescore, matchscore)

            # insertions: updated in place, don't depend on predecessors
            insertions(i, l2, scores, inserted)
            backGrphIdx[i+1, 1:] = numpy.where(inserted, i+1, backGrphIdx[i+1, 1:])
            backStrIdx[i+1, 1:] = numpy.where(inserted, numpy.arange(l2), backStrIdx[i+1, 1:])

            # if we're doing local alignment, don't let bad global alignment
            # drag us negative
            if not self.globalAlign:
                backGrphIdx[i+1, :] = numpy.where(scores[i+1, :] > 0, backGrphIdx[i+1, :], -1)
                backStrIdx [i+1, :] = numpy.where(scores[i+1, :] > 0, backStrIdx[i+1, :], -1)
                scores[i+1, :]      = numpy.maximum(scores[i+1, :], 0)

        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)

    def prevIndices(self, node, nodeIDtoIndex):
        """Return a list of the previous dynamic programming table indices
           corresponding to predecessors of the current node."""
        prev = [nodeIDtoIndex[predID] for predID in list(node.inEdges.keys())]
        # if no predecessors, point to just before the graph
        if not prev:
            prev = [-1]
        return prev

    def initializeDynamicProgrammingData(self):
        """Initalize the dynamic programming tables:
            - set up scores array
            - set up backtracking array
            - create index to Node ID table and vice versa"""
        l1 = self.graph.nNodes
        l2 = len(self.sequence)

        nodeIDtoIndex = {}
        nodeIndexToID = {-1: None}
        # generate a dict of (nodeID) -> (index into nodelist (and thus matrix))
        ni = self.graph.nodeiterator()
        for (index, node) in enumerate(ni()):
            nodeIDtoIndex[node.ID] = index
            nodeIndexToID[index] = node.ID

        # Dynamic Programming data structures; scores matrix and backtracking
        # matrix
        scores = numpy.zeros((l1+1, l2+1), dtype=int)

        # initialize insertion score
        # if global align, penalty for starting at head != 0
        if self.globalAlign:
            scores[0, :] = numpy.arange(l2+1)*self._gap

            ni = self.graph.nodeiterator()
            for (index, node) in enumerate(ni()):
                prevIdxs = self.prevIndices(node, nodeIDtoIndex)
                best = scores[prevIdxs[0]+1, 0]
                for prevIdx in prevIdxs:
                    best = max(best, scores[prevIdx+1, 0])
                scores[index+1, 0] = best + self._gap

        # backtracking matrices
        backStrIdx = numpy.zeros((l1+1, l2+1), dtype=int)
        backGrphIdx = numpy.zeros((l1+1, l2+1), dtype=int)

        return nodeIDtoIndex, nodeIndexToID, scores, backStrIdx, backGrphIdx

    def backtrack(self, scores, backStrIdx, backGrphIdx, nodeIndexToID):
        """Backtrack through the scores and backtrack arrays.
           Return a list of sequence indices and node IDs (not indices, which
           depend on ordering)."""
        besti, bestj = scores.shape
        besti -= 1
        bestj -= 1
        if not self.globalAlign:
            besti, bestj = numpy.argwhere(scores == numpy.amax(scores))[-1]
        else:
            ni = self.graph.nodeiterator()
            # still have to find best final index to start from
            terminalIndices = [
                index for (index, node) in enumerate(ni()) if node.outDegree == 0
            ]

            besti = terminalIndices[0] + 1
            bestscore = scores[besti, bestj]
            for i in terminalIndices[1:]:
                score = scores[i+1, bestj]
                if score > bestscore:
                    bestscore, besti = score, i+1

        matches = []
        strindexes = []
        while ((self.globalAlign or scores[besti, bestj] > 0)) and (
            besti != 0 or bestj != 0
        ):
            nexti, nextj = backGrphIdx[besti, bestj], backStrIdx[besti, bestj]
            curstridx, curnodeidx = bestj-1, nodeIndexToID[besti-1]

            strindexes.insert(0, curstridx if nextj != bestj else None)
            matches.insert   (0, curnodeidx if nexti != besti else None)

            besti, bestj = nexti, nextj

        return strindexes, matches
    
"""
Tests for Issue #6, provided by @rlorigro
"""

def generate_poa_graph(sequences):
    """
    Initialize graph and align all sequences
    :param sequences: sequences to align
    :return: graph: the completed POA graph resulting from the given sequences
    """
    init_sequence = sequences[0]
    init_label = "0"

    graph = POAGraph(init_sequence, init_label)

    for i in range(1, len(sequences)):
        sequence = sequences[i]
        alignment = SeqGraphAlignment(sequence, graph,
                                                        fastMethod=False,
                                                        globalAlign=True,
                                                        matchscore=1,
                                                        mismatchscore=-1,
                                                        gapscore=-2)

        graph.incorporateSeqAlignment(alignment, sequence, str(i))

    return graph


def sequences_and_test(sequences, test_sequence):
    graph = generate_poa_graph(sequences)
    alignment = SeqGraphAlignment(test_sequence, graph,
                                                    fastMethod=False,
                                                    globalAlign=True,
                                                    matchscore=1,
                                                    mismatchscore=-1,
                                                    gapscore=-2)

    graph.incorporateSeqAlignment(alignment, test_sequence, "test")
    alignments = graph.generateAlignmentStrings()

    result = alignments[-2][1].replace("-", "")
    return graph, result


def test_order_of_alignment_case1():
    sequences = ['ATATTGTGTAAGGCACAATTAACA',
                 'ATATTGCAAGGCACAATTCAACA',
                 'ATATTGCAAGGCACACAACA',
                 'ATGTGCAAGAGCACATAACA']
    test_sequence = "ATATTGCAAGGCACACTAACA"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence


def test_order_of_alignment_case2():
    sequences = ["TTA", "TGC"]
    test_sequence = "TTGC"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence


def test_order_of_alignment_case3():
    sequences = ["TAGTGAAAGAGGAAAAGAA"]
    test_sequence = "GCCCAGAAATTCCAGACCAGC"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence


def test_order_of_alignment_case4():
    sequences = ["CTACTTGGGAGGCTGAGGTGG", "CCACTTGAGTTGAGG", "CTACTTGGGAAGCTAGAGGTGG"]
    test_sequence = "CTACTTGGGAGGCTGAGGTGG"

    _, result = sequences_and_test(sequences, test_sequence)
    assert result == test_sequence