from unittest import TestCase
import unittest

from  blockassembly.graph.graph import graph_multi_k

import numpy as np


spades = {"ref_seq":'ggacatcagata',
          "reads":["acat", "catc", "atca", "tcag", "caga", "agat", "gata", "tagg", "ggac"]}
multik = {"ref_seq":"aaaatcgatctcatcgaatt",
          "reads":["aaaatcgatctc","tctcatcgaatt"]}
    
param_list = [
            (["exp_mini",{"ref_seq":"ctata", "reads":["ctat","tata"], "shuffle":False}, [4], [4]],
              [[[],
                ['ctata'],
                []]]),
            (["exp_spades",{"ref_seq":spades["ref_seq"], "reads":spades["reads"], "shuffle":True}, [3,4,3], [3,4,4]],
              [[[(0, 1), (0, 9), (1, 2), (1, 11), (2, 3), (2, 7), (3, 4), (3, 5), (4, 10), (5, 0), (6, 3), (6, 7), (7, 8), (8, 1), (8, 9), (9, 4), (9, 5), (10, 2), (10, 11), (11, 6)],
                ['atca', 'cag', 'aga', 'gat', 'atag', 'agga', 'gaca', 'cat'],
                [(0, 1), (0, 7), (1, 2), (1, 5), (2, 3), (2, 6), (3, 0), (3, 4), (4, 2), (4, 5), (5, 3), (5, 6), (6, 1), (6, 7), (7, 0), (7, 4)]],
                [[(0, 1), (1, 7), (3, 0), (5, 6), (6, 3), (7, 2)],
                ['acatcagata', 'ggac', 'tagg'],
                []],
                [[(0, 1), (1, 7), (2, 9), (3, 0), (4, 11), (5, 6), (6, 3), (7, 2), (8, 10), (9, 8), (10, 4), (11, 5)],
                ['cagataggacatcag'],
                [(0, 0)]]]),
            (["exp_multik",{"ref_seq":multik["ref_seq"], "reads":multik["reads"], "shuffle":True} , [5,7,5], [5,7,7]],
            [[[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 11), (5, 6), (6, 7), (8, 9), (9, 10), (10, 4), (11, 12), (12, 13), (13, 14), (14, 0)],
                ['tcgatctcatcg', 'atcga', 'tcgaatt', 'aaaatcg'],
                [(0, 1), (1, 0), (1, 2), (3, 1)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11)],
                ['tctcatcgaatt', 'aaaatcgatctc'],
                []],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 0)],
                ['aaaatcgatctcatcgaatt'],
                []]]),
            (["exp1",{"ref_seq":'cddcdbdcbabbbdaadcccdccddaccdbcdbcdbdabcdadd', "read_length":15, "mean_coverage":5} , [4,5,6,7,8,4], [4,5,6,7,8,8]],
            [[[(0, 1), (1, 2), (2, 3), (2, 28), (3, 4), (3, 10), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 29), (10, 11), (10, 31), (11, 12), (12, 13), (13, 14), (13, 33), (14, 15), (15, 16), (17, 18), (18, 19), (18, 22), (18, 27), (19, 20), (20, 17), (20, 21), (21, 19), (21, 22), (21, 27), (22, 0), (22, 23), (23, 24), (24, 25), (25, 26), (26, 19), (26, 22), (26, 27), (27, 3), (27, 28), (28, 34), (29, 30), (30, 11), (30, 31), (31, 32), (32, 35), (33, 3), (33, 28), (34, 14), (34, 33), (35, 36), (36, 17), (36, 21)],
                ['ccdd', 'cddcdb', 'cdbd', 'dbdcbabbbda', 'bdabcd', 'bcdadd', 'dcccd', 'ccdcc', 'dccd', 'cddaccd', 'ccdb', 'cdbcd', 'bdaadcc', 'bcdb'],
                [(0, 1), (0, 9), (1, 2), (1, 11), (2, 3), (3, 4), (3, 12), (4, 5), (4, 13), (6, 0), (6, 7), (6, 10), (7, 6), (7, 8), (8, 0), (8, 7), (8, 10), (9, 0), (9, 7), (9, 10), (10, 2), (10, 11), (11, 5), (11, 13), (12, 6), (12, 8), (13, 2), (13, 11)]],
                [[(0, 1), (1, 2), (2, 3), (2, 9), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 27), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (16, 17), (17, 18), (18, 19), (19, 20), (20, 21), (21, 22), (22, 23), (23, 24), (24, 25), (25, 26), (26, 32), (27, 28), (28, 29), (29, 30), (31, 3), (31, 9), (32, 33), (33, 31), (33, 34), (34, 32), (35, 36), (36, 16)],
                ['cddcdbd', 'cdbdcbabbbdaad', 'cdbdabcdadd', 'aadcccdccddaccdbc', 'cdbcdb'],
                [(0, 1), (0, 2), (3, 4)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 24), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20), (20, 21), (21, 22), (22, 23), (23, 29), (24, 25), (25, 26), (26, 27), (28, 8), (29, 30), (30, 31), (30, 33), (31, 32), (32, 30), (33, 28), (34, 35), (35, 14)],
                ['cddcdbdcbabbbdaad', 'dbcdbdabcdadd', 'cdbcdb', 'aadcccdccddaccdbcd', 'dbcdbcd'],
                [(2, 1), (3, 2), (2, 4), (4, 2)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 21), (7, 8), (8, 9), (9, 10), (10, 11), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20), (20, 26), (21, 22), (22, 23), (23, 24), (25, 7), (26, 27), (27, 28), (27, 31), (28, 29), (29, 30), (30, 28), (30, 31), (31, 32), (32, 25), (33, 34), (34, 12)],
                ['cddcdbdcbabbbdaad', 'cdbcdbdabcdadd', 'aadcccdccddaccdbcdb', 'cdbcdbcdb'],
                [(2, 1), (2, 3), (3, 1), (3, 3)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 18), (6, 7), (7, 8), (8, 9), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 23), (18, 19), (19, 20), (20, 21), (22, 6), (23, 24), (24, 25), (25, 26), (26, 27), (27, 28), (28, 29), (29, 30), (30, 22), (31, 32), (32, 10)],
                ['cddcdbdcbabbbdaad', 'aadcccdccddaccdbcdbcdbdabcdadd'],
                []],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 18), (6, 7), (7, 8), (8, 9), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 23), (18, 19), (19, 20), (20, 21), (21, 33), (22, 6), (23, 24), (24, 25), (25, 26), (26, 27), (27, 28), (28, 29), (29, 30), (30, 22), (31, 32), (32, 10), (33, 34), (34, 35), (35, 36), (36, 31)],
                ['cddcdbdcbabbbdaadcccdccddaccdbcdbcdbdabcdadd'],
                []]]),
            (["exp2",{"ref_seq":"jeabjabijaigedaegibiebdgfdjgjbjecghiijcaghibhbeaifehiicgciiggfgaagjbijbcijjfachd", "read_length":10, "mean_coverage":20} , [3,4,5,3], [3,4,5,5]],
            [[[(0, 1), (0, 51), (1, 2), (2, 3), (3, 4), (4, 5), (4, 25), (5, 6), (5, 36), (6, 7), (7, 31), (8, 0), (8, 9), (9, 10), (9, 15), (9, 45), (10, 11), (11, 12), (11, 43), (12, 13), (13, 14), (13, 23), (14, 10), (14, 15), (14, 45), (15, 37), (15, 39), (16, 17), (17, 18), (17, 21), (17, 64), (18, 19), (18, 52), (19, 20), (19, 28), (19, 34), (19, 46), (20, 18), (20, 21), (20, 64), (21, 22), (22, 14), (22, 23), (23, 20), (23, 28), (23, 34), (23, 46), (24, 5), (24, 25), (25, 26), (25, 33), (26, 27), (26, 65), (27, 32), (27, 35), (28, 29), (29, 30), (30, 58), (31, 8), (32, 26), (32, 33), (33, 19), (33, 52), (34, 32), (34, 35), (35, 6), (35, 36), (36, 37), (36, 39), (37, 38), (38, 70), (39, 40), (40, 41), (40, 57), (41, 59), (42, 12), (42, 43), (43, 44), (44, 0), (44, 9), (45, 20), (45, 28), (45, 34), (45, 46), (46, 47), (47, 48), (48, 16), (48, 49), (49, 44), (50, 1), (50, 51), (51, 19), (51, 52), (52, 53), (53, 54), (54, 55), (55, 56), (56, 41), (56, 57), (57, 67), (58, 68), (59, 60), (60, 61), (61, 16), (61, 49), (62, 63), (63, 17), (64, 27), (64, 65), (65, 24), (65, 66), (66, 42), (67, 62), (68, 69), (70, 71), (71, 72), (72, 73), (73, 50)],
                ['hib', 'ibhbea', 'eai', 'aifehi', 'hii', 'iicg', 'cgci', 'cii', 'iig', 'agj', 'gjb', 'jbi', 'bij', 'ijb', 'jbci', 'cij', 'bje', 'eab', 'abj', 'bja', 'ijjfachd', 'jab', 'abi', 'ija', 'jai', 'aig', 'igedaegib', 'iggf', 'gfgaag', 'jecg', 'cgh', 'ghi', 'iij', 'ijcag', 'agh', 'ibi', 'biebdgf', 'gfdjgj', 'jbj'],
                [(0, 1), (0, 35), (1, 2), (1, 17), (2, 3), (2, 25), (3, 0), (3, 4), (4, 5), (4, 8), (4, 32), (5, 6), (5, 30), (6, 7), (6, 15), (7, 5), (7, 8), (7, 32), (8, 26), (8, 27), (9, 10), (10, 11), (10, 14), (10, 38), (11, 12), (11, 36), (12, 13), (12, 20), (12, 23), (12, 33), (13, 11), (13, 14), (13, 38), (14, 7), (14, 15), (15, 13), (15, 20), (15, 23), (15, 33), (16, 29), (17, 18), (17, 22), (18, 16), (18, 19), (19, 21), (19, 24), (21, 18), (21, 22), (22, 12), (22, 36), (23, 21), (23, 24), (24, 3), (24, 25), (25, 26), (25, 27), (26, 1), (26, 35), (27, 28), (27, 37), (28, 9), (28, 34), (29, 6), (29, 30), (30, 31), (31, 0), (31, 4), (32, 13), (32, 20), (32, 23), (32, 33), (33, 9), (33, 34), (34, 31), (35, 12), (35, 36), (36, 28), (36, 37), (37, 10), (38, 16), (38, 19)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 27), (7, 8), (7, 42), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 36), (14, 15), (14, 63), (15, 16), (16, 17), (16, 31), (17, 18), (18, 19), (19, 20), (20, 24), (21, 22), (22, 23), (23, 49), (24, 25), (25, 26), (26, 57), (27, 28), (28, 7), (29, 30), (30, 17), (30, 31), (31, 32), (32, 33), (33, 34), (34, 35), (35, 76), (36, 37), (37, 38), (38, 75), (39, 40), (40, 41), (40, 47), (41, 8), (41, 42), (42, 43), (43, 44), (44, 45), (45, 48), (46, 41), (46, 47), (47, 0), (48, 46), (49, 29), (50, 51), (51, 52), (52, 53), (53, 54), (54, 55), (55, 56), (56, 67), (57, 69), (58, 59), (59, 60), (60, 14), (61, 62), (62, 15), (62, 63), (63, 64), (64, 65), (65, 66), (66, 39), (67, 68), (68, 61), (69, 70), (71, 72), (72, 73), (73, 74), (74, 50), (75, 58), (76, 71)],
                ['ghibhbeaifehii', 'hiicgciiggfgaagjb', 'gjbij', 'bijbcijjfachd', 'jeabjabij', 'bijaigedaegibiebdgfdjgjb', 'gjbjecghi', 'ghii', 'hiijcaghi'],
                [(0, 1), (0, 8), (1, 2), (1, 6), (2, 3), (2, 5), (4, 3), (4, 5), (5, 2), (5, 6), (6, 0), (6, 7), (7, 1), (7, 8), (8, 0), (8, 7)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 23), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 32), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 20), (18, 19), (19, 45), (20, 21), (21, 22), (22, 53), (23, 24), (24, 25), (25, 6), (26, 27), (27, 28), (28, 29), (29, 30), (30, 31), (31, 74), (32, 33), (33, 34), (34, 71), (35, 36), (36, 37), (37, 38), (38, 39), (39, 40), (40, 43), (41, 42), (42, 0), (43, 44), (44, 41), (45, 46), (46, 26), (47, 48), (48, 49), (49, 50), (50, 51), (51, 52), (52, 72), (53, 65), (54, 55), (55, 56), (56, 12), (57, 58), (58, 59), (59, 60), (60, 61), (61, 62), (62, 35), (63, 64), (64, 57), (65, 66), (67, 68), (68, 69), (69, 70), (70, 47), (71, 73), (72, 63), (73, 54), (74, 75), (75, 67)],
                ['jeabjabijaigedaegibiebdgfdjgjbjecghiijcaghibhbeaifehiicgciiggfgaagjbijbcijjfachd'],
                []],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 23), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 32), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 20), (18, 19), (19, 45), (20, 21), (21, 22), (22, 53), (23, 24), (24, 25), (25, 6), (26, 27), (27, 28), (28, 29), (29, 30), (30, 31), (31, 74), (32, 33), (33, 34), (34, 71), (35, 36), (36, 37), (37, 38), (38, 39), (39, 40), (40, 43), (41, 42), (42, 0), (43, 44), (44, 41), (45, 46), (46, 26), (47, 48), (48, 49), (49, 50), (50, 51), (51, 52), (52, 72), (53, 65), (54, 55), (55, 56), (56, 12), (57, 58), (58, 59), (59, 60), (60, 61), (61, 62), (62, 35), (63, 64), (64, 57), (65, 66), (67, 68), (68, 69), (69, 70), (70, 47), (71, 73), (72, 63), (73, 54), (74, 75), (75, 67)],
                ['jeabjabijaigedaegibiebdgfdjgjbjecghiijcaghibhbeaifehiicgciiggfgaagjbijbcijjfachd'],
                []]]),
            (["exp3",{"ref_seq":"abcdefghituvwxyzjkldefmnop", "read_length":10, "mean_coverage":20} , [3,4,5], [3,4,5]],
            [[[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 22), (8, 9), (9, 10), (9, 17), (10, 11), (11, 12), (12, 13), (14, 15), (15, 16), (16, 9), (17, 18), (18, 19), (19, 20), (20, 21), (21, 0), (22, 8)],
                ['efghituvwxyzjklde', 'def', 'efmnop', 'abcde'],
                [(0, 1), (1, 0), (1, 2), (3, 1)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 21), (7, 8), (7, 15), (8, 9), (9, 10), (10, 11), (12, 13), (13, 14), (14, 8), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20), (20, 0), (21, 22), (22, 7)],
                ['defghituvwxyzjkldef', 'defmnop', 'abcdef'],
                [(0, 0), (0, 1), (2, 0), (2, 1)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 19), (6, 7), (7, 8), (8, 9), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 0), (19, 20), (20, 21), (21, 6)],
                ['abcdefghituvwxyzjkldefmnop'],
                []]]),
            (["exp4",{"ref_seq":"abcdefghituvwxyzjkldefmntuvwxyzop", "read_length":10, "mean_coverage":20} , [3,4,5], [3,4,5]],
            [[[(0, 1), (1, 2), (2, 3), (3, 4), (3, 9), (4, 5), (5, 6), (6, 7), (7, 15), (8, 3), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 16), (15, 13), (16, 17), (17, 18), (18, 19), (18, 23), (19, 20), (21, 22), (22, 8), (23, 24), (24, 21)],
                ['abcde', 'def', 'efghitu', 'yzjklde', 'efmntu', 'tuvwxyz', 'yzop'],
                [(0, 1), (1, 2), (1, 4), (2, 5), (3, 1), (4, 5), (5, 3), (5, 6)]],
                [[(0, 1), (1, 2), (2, 3), (2, 8), (3, 4), (4, 5), (5, 6), (6, 14), (7, 3), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 16), (14, 15), (15, 13), (16, 17), (17, 18), (18, 19), (18, 23), (19, 20), (21, 22), (22, 7), (23, 24), (24, 25), (25, 21)],
                ['abcdef', 'defghituv', 'xyzjkldef', 'defmntuv', 'tuvwxyz', 'xyzop'],
                [(0, 1), (0, 3), (1, 4), (2, 1), (2, 3), (3, 4), (4, 2), (4, 5)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 12), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 15), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (17, 22), (18, 19), (20, 21), (21, 6), (22, 23), (23, 24), (24, 25), (25, 20)],
                ['abcdefghituvw', 'wxyzjkldefmntuvw', 'tuvwxyz', 'wxyzop'],
                [(0, 2), (1, 2), (2, 1), (2, 3)]]]),
            (["exp5",{"ref_seq":"jeabjabijaigedaegibibonsaiigeebrepeatdgfbjbonsaiabrepeatidjgjbigbonsaiejebjabicg", "read_length":10, "mean_coverage":20} , [3,4,5,6,7,8,9,3], [3,4,5,6,7,8,9,9]],
            [[[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 16), (4, 31), (4, 54), (5, 6), (6, 7), (6, 22), (6, 27), (7, 26), (8, 9), (9, 10), (10, 11), (10, 21), (11, 12), (11, 37), (12, 13), (13, 14), (14, 15), (15, 34), (16, 17), (17, 18), (18, 19), (18, 20), (19, 24), (19, 53), (20, 11), (20, 21), (21, 7), (21, 22), (21, 27), (22, 23), (22, 42), (23, 25), (23, 30), (24, 23), (24, 42), (25, 7), (25, 22), (25, 27), (26, 8), (27, 28), (27, 45), (27, 47), (27, 55), (28, 29), (29, 25), (29, 30), (30, 5), (30, 16), (30, 31), (30, 54), (31, 32), (31, 48), (32, 33), (32, 51), (33, 57), (34, 35), (35, 0), (35, 36), (36, 28), (36, 45), (36, 47), (36, 55), (37, 38), (38, 39), (39, 40), (40, 41), (41, 23), (41, 42), (42, 0), (42, 36), (43, 44), (43, 46), (44, 28), (44, 45), (44, 47), (44, 55), (45, 44), (45, 46), (46, 1), (47, 32), (47, 48), (48, 49), (49, 1), (50, 32), (50, 48), (51, 52), (52, 24), (52, 53), (53, 26), (54, 50), (55, 56), (57, 58), (58, 59), (59, 60), (60, 43)],
                ['jbo', 'bonsai', 'aiab', 'abr', 'brepea', 'eat', 'atidjgjb', 'aieje', 'eab', 'abj', 'bja', 'ebj', 'jab', 'abi', 'bija', 'jai', 'aig', 'ige', 'jbi', 'atdgfbj', 'bjb', 'gedaegib', 'ibi', 'bib', 'ibo', 'big', 'igbo', 'aiig', 'geeb', 'ebr', 'bicg'],
                [(0, 1), (1, 2), (1, 7), (1, 16), (1, 27), (2, 3), (2, 9), (2, 13), (3, 4), (4, 5), (4, 8), (5, 6), (5, 19), (6, 0), (6, 18), (8, 3), (8, 9), (8, 13), (9, 10), (9, 20), (10, 12), (10, 15), (11, 10), (11, 20), (12, 3), (12, 9), (12, 13), (13, 14), (13, 23), (13, 25), (13, 30), (14, 12), (14, 15), (15, 2), (15, 7), (15, 16), (15, 27), (16, 17), (16, 26), (17, 21), (17, 28), (18, 14), (18, 23), (18, 25), (18, 30), (19, 10), (19, 20), (20, 0), (20, 18), (21, 22), (21, 24), (22, 14), (22, 23), (22, 25), (22, 30), (23, 22), (23, 24), (24, 1), (25, 17), (25, 26), (26, 1), (27, 17), (27, 26), (28, 11), (28, 29), (29, 4)]],
                [[(0, 1), (1, 2), (2, 3), (3, 4), (3, 14), (3, 55), (4, 5), (5, 6), (6, 24), (7, 8), (8, 9), (9, 10), (9, 36), (10, 11), (11, 12), (12, 13), (13, 33), (14, 15), (15, 16), (16, 17), (17, 21), (18, 19), (19, 20), (20, 23), (21, 22), (22, 23), (23, 26), (24, 25), (25, 7), (26, 27), (26, 57), (27, 28), (28, 29), (29, 30), (30, 31), (31, 32), (31, 51), (32, 64), (33, 34), (34, 35), (35, 63), (36, 37), (37, 38), (38, 39), (39, 42), (40, 41), (41, 0), (42, 40), (43, 44), (44, 45), (45, 46), (46, 1), (47, 48), (48, 49), (49, 1), (50, 32), (50, 51), (51, 52), (52, 53), (53, 54), (54, 25), (55, 56), (56, 50), (57, 58), (59, 60), (60, 61), (61, 62), (62, 43), (63, 47), (64, 59)],
                ['eatdgfbjbon', 'bonsai', 'saiabre', 'brepeat', 'eatidjgjbigbon', 'saiejebja', 'jeabja', 'bjabi', 'abijaige', 'igedaegibibon', 'saiige', 'igeebre', 'abicg'],
                [(0, 1), (1, 2), (1, 5), (1, 10), (2, 3), (3, 0), (3, 4), (4, 1), (5, 7), (6, 7), (7, 8), (7, 12), (8, 9), (8, 11), (9, 1), (10, 9), (10, 11), (11, 3)]],
                [[(0, 1), (1, 2), (2, 3), (2, 12), (2, 64), (3, 4), (4, 5), (5, 21), (6, 7), (7, 8), (7, 33), (8, 9), (9, 10), (10, 11), (11, 30), (12, 13), (13, 14), (14, 15), (15, 18), (16, 17), (17, 41), (18, 19), (19, 20), (20, 42), (21, 22), (22, 23), (23, 6), (24, 25), (25, 26), (26, 27), (27, 28), (28, 29), (29, 66), (30, 31), (31, 32), (32, 63), (33, 34), (34, 35), (35, 36), (36, 39), (37, 38), (38, 0), (39, 40), (40, 37), (41, 42), (42, 24), (42, 57), (43, 44), (44, 45), (45, 46), (46, 1), (47, 48), (48, 49), (49, 1), (50, 51), (51, 52), (52, 53), (53, 54), (54, 23), (55, 56), (56, 50), (57, 58), (59, 60), (60, 61), (61, 62), (62, 43), (63, 65), (64, 55), (65, 47), (66, 67), (67, 59)],
                ['peatdgfbjbons', 'bonsai', 'nsaiabrep', 'brepeat', 'peatidjgjbigbons', 'nsaiejebjab', 'jeabjab', 'bjabi', 'jabijaigedaegibibons', 'nsaiigeebrep', 'jabicg'],
                [(0, 1), (1, 2), (1, 5), (1, 9), (2, 3), (3, 0), (3, 4), (4, 1), (5, 7), (6, 7), (7, 8), (7, 10), (8, 1), (9, 3)]],
                [[(0, 1), (1, 2), (1, 10), (1, 64), (2, 3), (3, 4), (4, 18), (5, 6), (5, 30), (6, 7), (7, 8), (8, 9), (9, 27), (10, 11), (11, 12), (12, 13), (13, 15), (14, 38), (15, 16), (16, 17), (17, 44), (18, 19), (19, 20), (20, 21), (21, 5), (22, 23), (23, 24), (24, 25), (25, 26), (26, 70), (27, 28), (28, 29), (29, 62), (30, 31), (31, 32), (32, 33), (33, 36), (34, 35), (35, 0), (36, 37), (37, 54), (38, 39), (39, 53), (39, 63), (40, 41), (41, 42), (42, 43), (43, 1), (44, 53), (44, 63), (45, 46), (46, 47), (47, 1), (48, 49), (49, 50), (50, 51), (51, 52), (52, 21), (53, 22), (54, 34), (55, 56), (56, 48), (58, 59), (59, 60), (60, 61), (61, 40), (62, 66), (63, 57), (64, 65), (65, 55), (66, 67), (67, 45), (68, 69), (69, 58), (70, 68)],
                ['epeatdgfbjbonsa', 'bonsai', 'onsaiabrepe', 'brepeat', 'epeatidjgjbigbonsa', 'onsaiejebjabi', 'jeabjabi', 'bjabijaigedaegibibonsa', 'onsaiigeebrepe', 'bjabicg'],
                [(0, 1), (1, 2), (1, 5), (1, 8), (2, 3), (3, 0), (3, 4), (4, 1), (5, 7), (5, 9), (6, 7), (6, 9), (7, 1), (8, 3)]],
                [[(0, 1), (0, 8), (0, 62), (1, 2), (2, 3), (3, 15), (4, 5), (5, 6), (6, 7), (7, 23), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 41), (15, 16), (16, 17), (17, 18), (18, 34), (19, 20), (20, 21), (21, 22), (22, 71), (23, 24), (24, 25), (25, 59), (26, 27), (27, 28), (28, 29), (29, 32), (30, 31), (31, 0), (32, 33), (33, 51), (34, 4), (34, 26), (35, 36), (36, 67), (37, 38), (38, 39), (39, 40), (40, 1), (40, 8), (40, 62), (41, 60), (42, 43), (43, 44), (44, 1), (44, 8), (44, 62), (45, 46), (46, 47), (47, 48), (48, 50), (49, 19), (50, 34), (51, 52), (52, 30), (53, 54), (54, 45), (55, 56), (56, 57), (57, 58), (58, 37), (59, 64), (60, 61), (62, 63), (63, 68), (64, 65), (65, 66), (66, 42), (67, 49), (68, 53), (69, 70), (70, 55), (71, 72), (72, 69)],
                ['repeatdgfbjbonsai', 'bonsaiabrepea', 'brepeat', 'repeatidjgjbigbonsai', 'bonsaiejebjabicg', 'jeabjabijaigedaegibibonsai', 'bonsaiigeebrepea'],
                [(0, 1), (0, 4), (0, 6), (1, 2), (2, 0), (2, 3), (3, 1), (3, 4), (3, 6), (5, 1), (5, 4), (5, 6), (6, 2)]],
                [[(0, 1), (1, 2), (2, 42), (3, 4), (4, 5), (5, 18), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 34), (12, 13), (13, 14), (14, 28), (15, 16), (16, 17), (17, 68), (18, 19), (19, 20), (20, 54), (21, 22), (22, 23), (23, 26), (24, 25), (25, 0), (26, 27), (27, 45), (28, 29), (28, 48), (29, 3), (30, 62), (31, 32), (32, 33), (33, 57), (34, 55), (35, 36), (36, 37), (37, 6), (38, 39), (39, 40), (40, 43), (41, 15), (42, 12), (43, 44), (44, 29), (44, 48), (45, 46), (46, 47), (47, 24), (48, 21), (49, 50), (50, 38), (51, 52), (52, 53), (53, 71), (54, 59), (55, 56), (57, 58), (58, 72), (59, 60), (60, 61), (61, 70), (62, 69), (63, 49), (64, 65), (65, 51), (66, 67), (67, 64), (68, 66), (69, 41), (70, 35), (71, 31), (72, 63)],
                ['brepeatdgfbjbonsaiabrepeat', 'brepeatidjgjbigbonsaiejebjabicg', 'jeabjabijaigedaegibibonsaiigeebrepeat'],
                [(0, 0), (0, 1), (2, 0), (2, 1)]],
                [[(0, 1), (1, 44), (2, 3), (3, 14), (4, 5), (5, 52), (6, 7), (7, 25), (8, 9), (9, 21), (10, 11), (11, 61), (12, 13), (13, 43), (14, 12), (15, 16), (16, 19), (17, 18), (18, 0), (19, 20), (20, 34), (21, 22), (22, 53), (23, 24), (24, 46), (25, 51), (26, 27), (27, 48), (28, 29), (29, 32), (30, 10), (31, 8), (32, 33), (33, 38), (34, 35), (35, 36), (36, 37), (37, 17), (38, 65), (39, 40), (40, 28), (41, 42), (42, 70), (43, 56), (44, 31), (46, 47), (47, 67), (48, 4), (49, 50), (50, 63), (51, 45), (52, 6), (53, 2), (54, 69), (55, 39), (56, 49), (57, 58), (58, 41), (59, 60), (60, 57), (62, 30), (63, 64), (64, 26), (65, 15), (66, 23), (67, 68), (68, 55), (69, 62), (70, 66)],
                ['jaigedaegibibonsaiigeebrepeatdgfbjbonsaiabrepeatidjgjbigbonsaiejebjabicg', 'jeabjabijaigeda'],
                []],
                [[(0, 1), (1, 44), (2, 3), (3, 14), (4, 5), (5, 52), (6, 7), (7, 25), (8, 9), (9, 21), (10, 11), (11, 61), (12, 13), (13, 43), (14, 12), (15, 16), (16, 19), (17, 18), (18, 0), (19, 20), (20, 34), (21, 22), (22, 53), (23, 24), (24, 46), (25, 51), (26, 27), (27, 48), (28, 29), (29, 32), (30, 10), (31, 8), (32, 33), (33, 38), (34, 35), (35, 36), (36, 37), (37, 17), (38, 65), (39, 40), (40, 28), (41, 42), (42, 70), (43, 56), (44, 31), (46, 47), (47, 67), (48, 4), (49, 50), (50, 63), (51, 45), (52, 6), (53, 2), (54, 69), (55, 39), (56, 49), (57, 58), (58, 41), (59, 60), (60, 57), (61, 71), (62, 30), (63, 64), (64, 26), (65, 15), (66, 23), (67, 68), (68, 55), (69, 62), (70, 66), (71, 59)],
                ['jeabjabijaigedaegibibonsaiigeebrepeatdgfbjbonsaiabrepeatidjgjbigbonsaiejebjabicg'],
                []]]),
            ]

class TestMultikAssemblyOnDifferentCases(TestCase):
    ## TODO add test on kmers and DBG edges not only on unitigs and compacted DBG
    def test_works_as_expected(self):
        for (name, kwargs,kmins,kmaxs), res in param_list:
            print("Testing {} ...".format(name))
            for k, (kmin, kmax) in enumerate(zip(kmins, kmaxs)):
                lg, lunitigs, lc_g = res[k]
                _, _, _, g, unitigs, c_g = graph_multi_k(**kwargs, kmin = kmin, kmax = kmax, verbose=False)
                idu = np.argsort(unitigs)
                unitigs = list(np.array(unitigs)[idu])
                c_e = list((tuple(int(n) for n in e)) for e in c_g.edges())
                c_e = [(list(idu).index(e[0]),list(idu).index(e[1])) for e in c_e]
                c_e.sort()
                lidu = np.argsort(lunitigs)
                lunitigs = list(np.array(lunitigs)[lidu])
                lc_e = [(list(lidu).index(e[0]),list(lidu).index(e[1])) for e in lc_g]
                lc_e.sort()
                
                for p1,p2 in zip([unitigs, c_e],[lunitigs, lc_e]):
                    with self.subTest((p1, p2)):
                        self.assertEqual(p1, p2)
                
if __name__=="__main__":
    unittest.main()