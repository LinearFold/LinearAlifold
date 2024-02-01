#!/usr/bin/env python2

from __future__ import division
from hopcroftkarp import HopcroftKarp
from collections import defaultdict
import sys


def pairs(s):
    stack1, stack2, stack3, stack4 = [], [], [], []
    ret = set()
    for i, x in enumerate(s, 1):
        if x == '(':
            stack1.append(i)
        elif x == ')':
            ret.add((stack1.pop(), i))
        elif x == '[':
            stack2.append(i)
        elif x == ']':
            ret.add((stack2.pop(), i))
        elif x == '{':
            stack3.append(i)
        elif x == '}':
            ret.add((stack3.pop(), i))
        elif x == '<':
            stack4.append(i)
        elif x == '>':
            ret.add((stack4.pop(), i))
    return ret

def eval(gold_pairs, test_pairs):
    # if len(gold_pairs) == 0:
    #     if len(test_pairs) == 0:
    #         return 1.0, 1.0, 1.0, 1.0
    #     else:
    #         return 0.0, 1.0, 0.0, 1.0
    slip_matched, non_slip_matched = 0, 0
    gold, test = len(gold_pairs), len(test_pairs)
    # print gold
    # print test
    graph = defaultdict(set)
    # graph = dict()
    for (i, j) in test_pairs:
        # non slip
        if (i,j) in gold_pairs:
            # print i,j
            non_slip_matched += 1
        # # slip
        # for (x,y) in [(i, j), (i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
        #     if (x,y) in gold_pairs:
        #         slip_matched += 1
        #         gold_pairs.remove((x,y))
        #         break

    ###########
        # # slip
        # graph[(i,j)] = set()
        for (x,y) in [(i, j), (i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
        # for (x,y) in [(i, j)]:
            if (x, y) in gold_pairs:
                # if (i,j) not in graph:
                #     graph[(i,j)] = set()
                graph[(i,j)].add((str(x), str(y)))

    # print(graph)

    # print HopcroftKarp(graph).maximum_matching()

    slip_matched = len(HopcroftKarp(graph).maximum_matching()) / 2

    ###########


    return slip_matched / (test + 1e-6), slip_matched / (gold + 1e-6), non_slip_matched / (test + 1e-6), non_slip_matched / (gold + 1e-6)


if __name__ == '__main__':

    golds, tests = open(sys.argv[1]).readlines(), open(sys.argv[2]).readlines()
    if len(golds) != len(tests):
        print('Gold file and test file are not matched!!!')
        exit()
    # tot = FScore()
    tot_P_slip, tot_R_slip, tot_P_non_slip, tot_R_non_slip = 0, 0, 0, 0
    for i, (gold_line, test_line) in enumerate(zip(golds, tests)):
        if len(gold_line) != len(test_line):
            print('gold %s' % (gold_line))
            print('test %s' % (test_line))
            print('gold length and test length are not matched!!!')
            exit()
        P_slip, R_slip, P_non_slip, R_non_slip = eval(pairs(gold_line), pairs(test_line))

        print('(P_slip = {:0.2f}, R_slip = {:0.2f}, P_noslip = {:0.2f}, R_noslip = {:0.2f})'.format(
            P_slip * 100,
            R_slip * 100,
            P_non_slip * 100,
            R_non_slip * 100))

        tot_P_slip += P_slip
        tot_R_slip += R_slip
        tot_P_non_slip += P_non_slip
        tot_R_non_slip += R_non_slip

    print('total (P_slip = {:0.2f}, R_slip = {:0.2f}, P_noslip = {:0.2f}, R_noslip = {:0.2f})'.format(
            tot_P_slip * 100 / len(golds),
            tot_R_slip * 100 / len(golds),
            tot_P_non_slip * 100 / len(golds),
            tot_R_non_slip * 100 / len(golds)))
    