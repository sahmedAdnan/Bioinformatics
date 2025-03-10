import unittest
import numpy as np

from linear.alignment import global_linear_alignment
from affine.alignment import global_affine_alignment

class TestProject2Examples(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Define the scoring matrix and store the sequences/expected values
        from project2_examples.txt.
        """
        cls.score_matrix = {
            'A': {'A': 0, 'C': 5, 'G': 2, 'T': 5},
            'C': {'A': 5, 'C': 0, 'G': 5, 'T': 2},
            'G': {'A': 2, 'C': 5, 'G': 0, 'T': 5},
            'T': {'A': 5, 'C': 2, 'G': 5, 'T': 0}
        }

        # We convert sequences to uppercase so we can index properly into the above matrix
        # Case 1:
        cls.seq1_case1 = "acgtgtcaacgt".upper() 
        cls.seq2_case1 = "acgtcgtagcta".upper()  

        cls.linear_score_case1 = 22
        cls.affine_score_case1 = 24
        cls.num_alignments_linear_case1 = 2
        

        # Case 2:
        cls.seq1_case2 = "aataat".upper()  
        cls.seq2_case2 = "aagg".upper()    

        cls.linear_score_case2 = 14
        cls.affine_score_case2 = 22
        cls.num_alignments_linear_case2 = 1  
        cls.num_alignments_affine_case2 = 3  

        # Case 3:
        cls.seq1_case3 = "tccagaga".upper()  
        cls.seq2_case3 = "tcgat".upper()     

        cls.linear_score_case3 = 20
        cls.affine_score_case3 = 29
        cls.num_alignments_linear_case3 = 4
        cls.num_alignments_affine_case3 = 1

        # Case 4: large sequences, we only verify the final score
        cls.seq1_case4 = (
            "ggcctaaaggcgccggtctttcgtaccccaaaatctcggcattttaagataagtgagtgttgcgttacactagcgatcta"
            "ccgcgtcttatacttaagcgtatgcccagatctgactaatcgtgcccccggattagacgggcttgatgggaaagaacagc"
            "tcgtctgtttacgtataaacagaatcgcctgggttcgc"
        ).upper()
        cls.seq2_case4 = (
            "gggctaaaggttagggtctttcacactaaagagtggtgcgtatcgtggctaatgtaccgcttctggtatcgtggcttacg"
            "gccagacctacaagtactagacctgagaactaatcttgtcgagccttccattgagggtaatgggagagaacatcgagtca"
            "gaagttattcttgtttacgtagaatcgcctgggtccgc"
        ).upper()

        cls.linear_score_case4 = 325
        cls.affine_score_case4 = 395

        # Gap penalties for linear: g(k)=5*k => gap_penalty=5
        cls.linear_gap_penalty = 5
        # Gap penalties for affine: g(k)=5 + 5*k => gap_open=5, gap_extend=5
        cls.gap_open = 5
        cls.gap_extend = 5

    # ------------------------------------------------------
    # Case 1: seq1=ACGTGTCAACGT, seq2=ACGTCGTAGCTA
    # ------------------------------------------------------
    def test_case1_linear(self):
        score, alignments = global_linear_alignment(
            self.seq1_case1,
            self.seq2_case1,
            self.score_matrix,
            gap_penalty=self.linear_gap_penalty,
            backtracking=True
        )
        self.assertEqual(score, self.linear_score_case1,
                         "Case1 (linear): score mismatch.")
        self.assertEqual(len(alignments), self.num_alignments_linear_case1,
                         "Case1 (linear): number of optimal alignments mismatch.")

    def test_case1_affine(self):
        score, alignments = global_affine_alignment(
            self.seq1_case1,
            self.seq2_case1,
            self.score_matrix,
            gap_open=self.gap_open,
            gap_extend=self.gap_extend,
            backtracking=True,
            show_all=True
        )
        self.assertEqual(score, self.affine_score_case1,
                         "Case1 (affine): score mismatch.")
        self.assertEqual(len(alignments), 1,
                         "Case1 (affine): number of alignments should be 1.")

    # ------------------------------------------------------
    # Case 2: seq1=AATAAT, seq2=AAGG
    # ------------------------------------------------------
    def test_case2_linear(self):
        score, alignments = global_linear_alignment(
            self.seq1_case2,
            self.seq2_case2,
            self.score_matrix,
            gap_penalty=self.linear_gap_penalty,
            backtracking=True
        )
        self.assertEqual(score, self.linear_score_case2,
                         "Case2 (linear): score mismatch.")
        self.assertEqual(len(alignments), self.num_alignments_linear_case2,
                         "Case2 (linear): number of alignments mismatch.")

    def test_case2_affine(self):
        score, alignments = global_affine_alignment(
            self.seq1_case2,
            self.seq2_case2,
            self.score_matrix,
            gap_open=self.gap_open,
            gap_extend=self.gap_extend,
            backtracking=True,
            show_all=True
        )
        self.assertEqual(score, self.affine_score_case2,
                         "Case2 (affine): score mismatch.")
        self.assertEqual(len(alignments), self.num_alignments_affine_case2,
                         "Case2 (affine): number of alignments mismatch.")

    # ------------------------------------------------------
    # Case 3: seq1=TCCAGAGA, seq2=TCGAT
    # ------------------------------------------------------
    def test_case3_linear(self):
        score, alignments = global_linear_alignment(
            self.seq1_case3,
            self.seq2_case3,
            self.score_matrix,
            gap_penalty=self.linear_gap_penalty,
            backtracking=True
        )
        self.assertEqual(score, self.linear_score_case3,
                         "Case3 (linear): score mismatch.")
        self.assertEqual(len(alignments), self.num_alignments_linear_case3,
                         "Case3 (linear): number of alignments mismatch.")

    def test_case3_affine(self):
        score, alignments = global_affine_alignment(
            self.seq1_case3,
            self.seq2_case3,
            self.score_matrix,
            gap_open=self.gap_open,
            gap_extend=self.gap_extend,
            backtracking=True,
            show_all=True
        )
        self.assertEqual(score, self.affine_score_case3,
                         "Case3 (affine): score mismatch.")
        self.assertEqual(len(alignments), self.num_alignments_affine_case3,
                         "Case3 (affine): number of alignments mismatch.")

    # ------------------------------------------------------
    # Case 4: large sequences
    # ------------------------------------------------------
    def test_case4_linear(self):
        score = global_linear_alignment(
            self.seq1_case4,
            self.seq2_case4,
            self.score_matrix,
            gap_penalty=self.linear_gap_penalty,
            backtracking=False
        )
        self.assertEqual(score, self.linear_score_case4,
                         "Case4 (linear): final score mismatch.")

    def test_case4_affine(self):
        score = global_affine_alignment(
            self.seq1_case4,
            self.seq2_case4,
            self.score_matrix,
            gap_open=self.gap_open,
            gap_extend=self.gap_extend,
            backtracking=False,
            show_all=False
        )
        self.assertEqual(score, self.affine_score_case4,
                         "Case4 (affine): final score mismatch.")


if __name__ == '__main__':
    unittest.main()
