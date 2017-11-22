#!/usr/bin/env python

from __future__ import division
from collections import defaultdict
import sys
from math import log
from itertools import product
startsym, stopsym = "<s>", "</s>"

def readfile(filename):
    for line in open(filename):
        wordtags = map(lambda x: x.rsplit("/", 1), line.split())
        yield [w for w,t in wordtags], [t for w,t in wordtags] # (word_seq, tag_seq) pair
    
def mle(filename): # Max Likelihood Estimation of HMM
    twfreq = defaultdict(lambda : defaultdict(int))
    ttfreq = defaultdict(lambda : defaultdict(int)) 
    tagfreq = defaultdict(int)    
    dictionary = defaultdict(set)

    for words, tags in readfile(filename):
        last = startsym
        tagfreq[last] += 1
        for word, tag in zip(words, tags) + [(stopsym, stopsym)]:
            #if tag == "VBP": tag = "VB" # +1 smoothing
            twfreq[tag][word] += 1            
            ttfreq[last][tag] += 1
            dictionary[word].add(tag)
            tagfreq[tag] += 1
            last = tag            
    
    model = defaultdict(float)
    num_tags = len(tagfreq)
    for tag, freq in tagfreq.iteritems(): 
        logfreq = log(freq)
        for word, f in twfreq[tag].iteritems():
            model[tag, word] = log(f) - logfreq 
        logfreq2 = log(freq + num_tags)
        for t in tagfreq: # all tags
            model[tag, t] = log(ttfreq[tag][t] + 1) - logfreq2 # +1 smoothing
        
    return dictionary, model

def liangs_decode(words, dictionary, model, gram=2):
    def backtrack(i, tag):
        if i == 0:
            return []
        return backtrack(i-1, back[i][tag]) + [tag]

    words = [startsym] + words + [stopsym]
    best = defaultdict(lambda: defaultdict(lambda: float("-inf")))
    best[0][startsym] = 1
    back = defaultdict(dict)

    for i, word in enumerate(words[1:], 1):
        for tag in dictionary[word]:
            if gram == 2:
                for prev in best[i-1]:
                    score = best[i-1][prev] + model[prev, tag] + model[tag, word]
                    if score > best[i][tag]:
                        best[i][tag] = score
                        back[i][tag] = prev
            elif gram == 3:
                if i == 1:
                    score = best[0][startsym] + model[startsym, tag] + model[tag, word]
                    back[1][tag] = '<s>'
                    if score > best[i][tag]:
                        best[i][startsym, tag] = score
                        back[i][tag] = startsym
                else:
                    print 'reading from best'
                    for prev in best[i-1]:
                        for prev_prev in best[i-2]:
                            score = best[i-1][prev_prev, prev] + model[prev_prev, prev, tag] + model[tag, word]
                        if score > best[i][tag]:
                            print 'writing to best'
                            best[i][prev, tag] = score
                            back[i][tag] = prev

    mytags = backtrack(len(words)-1, stopsym)[:-1]
    #print " ".join("%s/%s" % wordtag for wordtag in mywordtags)
    return mytags

def decode(words, dictionary, model, gram=2):
    words = [startsym] + words + [stopsym]
    dictionary[startsym] = [startsym]
    dictionary[stopsym] = [stopsym]
    n = len(words)
    pi = [defaultdict(lambda: float("-inf")) for w in words]
    pi[0][tuple(startsym for g in range(gram - 1))] = 1.

    final_cand = [tuple([t for t in dictionary[w]]) for w in words[n - gram + 1:]]
    zseq = list(max(product(*final_cand), key=lambda s : model[s] + sum(model[s[i], words[i]] for i in range(gram - 1))))

    def get_pi( j, subseq ):
        if subseq in pi[j]:
            return pi[j][subseq]
        elif j < gram - 1:
            return float( '-inf' )
        else:
            return max( [get_pi( j-1, (t,) + subseq[:-1] ) + model[(t,) + subseq] + model[words[j], subseq[-1]]
                         for t in dictionary[words[j-gram+1]]] )

    def backtrace(j, subseq):
        # return max(dictionary[words[j-gram+1]],
        #            key=lambda t : get_pi(j-1, (t,) + subseq[:-1]) + model[(t,) + subseq] + model[words[j], subseq[-1]])
        return max(dictionary[words[j-gram+1]],
                   key=lambda t : pi[j-1][(t,) + subseq[:-1]] + model[(t,) + subseq] + model[words[j], subseq[-1]])

    for i in range(n - 1, gram - 1, -1):
        # zseq.insert(0, backtrace(i, (zseq[-2], zseq[-1])))
        zseq.insert(0, backtrace(i, tuple(zseq[-gram+1:])))

    if startsym in zseq:
        zseq.remove(startsym)
    if stopsym in zseq:
        zseq.remove(stopsym)
    return zseq



def decode_tri(words, dictionary, model, gram=3):
    pi = defaultdict(lambda : float("-inf"))
    bp = {} #defaultdict( lambda: list() )
    words = [startsym] + words + [stopsym]
    pi[0, startsym, startsym] = 1.
    for t in dictionary[words[1]]:
        pi[1, startsym, t] = pi[0, startsym, startsym] + model[startsym, t] + model[t, words[1]]
    # bp[startsym, startsym] = []
    for i, word in enumerate(words[2:], 2):
        for prev in dictionary[words[i-1]]:
            for cur in dictionary[word]:
                p_prev = bp[i, prev, cur] = max(
                    dictionary[words[i - 2]],
                    key=lambda pp : pi[i-1, pp, prev] + model[pp, prev, cur] + model[cur, word]
                )
                pi[i, prev, cur] = pi[i-1, p_prev, prev] + model[p_prev, prev, cur] + model[cur, word]

    last = max(dictionary[words[-2]], key=lambda p : pi[len(words)-1, p, stopsym])

    def backtrack(i, p_prev, prev):
        if i == 1:
            return [prev]
        return backtrack(i-1, bp[i, p_prev, prev], p_prev) + [prev]

    return backtrack(len(words)-1, last, stopsym)[:-1]



def test(filename, dictionary, model, gram=2):
    errors = tot = 0
    for words, tags in readfile(filename):
        # mytags = liangs_decode(words, dictionary, model, gram)
        mytags = decode_tri(words, dictionary, model, gram)
        errors += sum(t1!=t2 for (t1,t2) in zip(tags, mytags))
        tot += len(words) 
        
    return errors/tot
        
if __name__ == "__main__":
    trainfile, devfile = sys.argv[1:3]
    dictionary, model = mle(trainfile)
    dictionary[startsym].add(startsym)
    dictionary[stopsym].add(stopsym)

    # print "train_err {0:.2%}".format(test(trainfile, dictionary, model))
    # print "dev_err {0:.2%}".format(test(devfile, dictionary, model))

    for words, tags in readfile(devfile):
        print words
        print decode_tri(words, dictionary, model, 3)
        print liangs_decode(words, dictionary, model)

    # sentences_labels = [wt for wt in readfile('dev.txt.lower.unk')]
    # for s, l in sentences_labels:
    #     for w in s:
    #         print w, '\t',
    #     print
    #     for t in l:
    #         print t, '\t',
    #     print
    #     my_labels = decoder(s, dictionary, model)
    #     liang = decode(s, dictionary,model)
    #     for i, ml in enumerate(my_labels):
    #         if ml != l[i]:
    #             print ml,'\t',
    #         else:
    #             print ' - \t',
    #     print
    #     for i, ll in enumerate(liang):
    #         if ll != l[i]:
    #             print ll, '\t',
    #         else:
    #             print ' - \t',
    #     print '\n'
    #
    #
    # print "dev_err {0:.2%}".format(test(devfile, dictionary, model))
