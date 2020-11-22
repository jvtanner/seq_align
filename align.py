"""
This file provides the skeleton code for align.py.

Locations with "FILL IN" in comments are where you need to add code. Implementational notes and hints
are incorporated throughout. DO NOT modify the static types that are pre-defined for each function unless
explicitly mentioned. Changing expected behavior of functions will likely result in loss of points by the
autograder.

Usage: python align.py input_file output_file
"""

import sys
from typing import Set, Tuple  # NOTE: You may need to "pip install typing" locally if this import gives you errors


#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a: float, b: float) -> bool:
    """
    Checks if two floating point numbers are equivalent.
    :param a: first number
    :param b: second number
    :return: True if a and be are within epsilon = 10^-6
    """
    epsilon = 10 ** (-6)
    return (abs(a - b) < epsilon)


#### ------- CLASSES ------- ####
class MatchMatrix(object):
    """
    A class representation of the match matrix, S. Stores the scores of matches between characters.
    """

    def __init__(self):
        """
        Initialize MatchMatrix class.

        NOTE: There is a bit of freedom in how to implement this.
        """

        self.matrix = {}

    def set_score(self, a: str, b: str, score: float) -> None:
        """
        Updates or adds a score for a specified match

        :param a: the character from sequence A
        :param b: the character from sequence B
        :param score: the score to set the match M(a,b)
        """
        ### FILL IN ###
        # Create an array using the input from line 9 till line (line5*line7)
        # should be size line5 * line7
        # Going left to right starting at line9: (i), (j), (seqA[i]), (seqB[j]), (match value)

        # Q: Doesn't make sense to use a and b as parameters--where are they coming from?
        # A: This is a method of MatchMatrix class. Another function can call this method,
        # and pass in parameter values for a and b.

        if (a, b) not in self.matrix:
            self.matrix[(a, b)] = score
        else:
            pass

    def get_score(self, a: str, b: str) -> float:
        """
        Returns the score for a particular match, where a is the
        character from sequence A and b is from sequence B.

        :param a: the character from sequence A
        :param b: the character from sequence B
        :return: the score of that match, M(a,b)
        """
        ### FILL IN ###
        return self.matrix[(a, b)]


class ScoreEntry(object):
    """
    A class object representing the score for each of the entries into one of
    the score matrices.
    """

    def __init__(self, ):
        self.score = 0
        self.pointers = set()


class ScoreMatrix(object):
    """
    A class representation of the score matrices (M, Ix, Iy), which will be dynamically updated.
    The score matrix consists of a 2-D array of ScoreEntries that are updated during alignment
    and used to output the maximum alignment.
    """

    def __init__(self, name: str, nrow: int, ncol: int) -> None:
        """
        Initialize ScoreMatrix class.

        :param name: identifier for the score matrix, should be in {Ix, Iy, M}
        :param nrow: number of rows for ScoreMatrix
        :param ncol: number of columns for ScoreMatrix
        """
        self.name = name
        self.nrow = nrow
        self.ncol = ncol

        self.array = []

        # Create an instance of the ScoreEntry class for each spot in the matrix
        for _ in range(nrow):
            row = [ScoreEntry() for _ in range(ncol)]
            self.array.append(row)

        # you need to figure out a way to represent this and how to initialize
        # Hint: it may be helpful to create a ScoreEntry class so that each entry in the
        # ScoreMatrix is an instance of this ScoreEntry class. This will help you keep track of all the
        # information that must be stored in each entry of the score matrix

    def get_score(self, row: int, col: int) -> float:
        """
        Return the current score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :return: current score at position (row, col)
        """

        ### FILL IN ###
        return self.array[row][col].score

    def set_score(self, row: int, col: int, score: float) -> None:
        """
        Set the score for an entry in ScoreMatrix.

        :param row: row index
        :param col: column index
        :param score: score to set at position (row, col)
        """

        ### FILL IN ###
        self.array[row][col].score = score

    def get_pointers(self, row: int, col: int) -> Set[Tuple[int, int, str]]:
        """
        Return the indices of the entries being pointed to. Remember, these are essential for the final traceback.

        :param row: row index
        :param col: column index
        :return: a set of indices (represented as tuples) corresponding to other entries being pointed to for traceback
        ex: {(1,0), (1,1)}
        """

        ### FILL IN ###
        return self.array[row][col].pointers

    def set_pointers(self, row: int, col: int, pointers: Set[Tuple[int, int, str]]) -> None:
        """
        Add pointers to each entry in the score matrix.

        :param row: row index
        :param col: column index
        :param pointers: set of pointers to add to your entry
        """

        ### FILL IN ###
        self.array[row][col].pointers = pointers

    def print_scores(self) -> str:
        """
        Returns a nicely formatted string containing the scores in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        :return: a string representation of the scores in the score matrix
        """

        ### FILL IN ###
        for i in range(self.nrow):
            for j in range(self.ncol):
                print(str(round(self.get_score(i, j), 2)) + ' ', end='')
            print('\n')

    def print_pointers(self) -> str:
        """
        Returns a nicely formatted string containing the pointers in the score matrix. This function is OPTIONAL
        (i.e. will not be checked by autograder) but will be extremely helpful for debugging!
        """

        pass


class AlignmentParameters(object):
    """
    A class to hold the alignment parameters.
    """

    def __init__(self) -> None:
        """
        Initialize AlignmentParameters object with default parameters.
        """

        # The definition for all of these class variables are documented in the annotated version of the input file
        # on the P1 page on Canvas.
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.len_alphabet_a = 0
        self.alphabet_a = ""
        self.len_alphabet_b = 0
        self.alphabet_b = ""
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file: str) -> None:
        """
        Read the parameters from an input file and update the alignment parameters accordingly

        :param input_file: path to the alignment input file (whose structure is defined on the project page)
        """

        ### FILL IN ###
        f = open(input_file, 'r')
        lines = [line.strip() for line in f.readlines()]

        # sequences
        self.seq_a = lines[0]
        self.seq_b = lines[1]

        # Global vs. Local
        if not int(lines[2]):
            self.global_alignment = True

        # Various penalties
        penalties = lines[3].split()
        self.dx = float(penalties[0])
        self.ex = float(penalties[1])
        self.dy = float(penalties[2])
        self.ey = float(penalties[3])

        self.len_alphabet_a = int(lines[4])
        self.len_alphabet_b = int(lines[6])

        self.alphabet_a = lines[5]
        self.alphabet_b = lines[7]

        # Parse the coordinates and associated values from the input file, then set scores in match matrix
        for line in lines[8:]:
            if line == '':
                continue
            _, _, a, b, score = line.split()
            self.match_matrix.set_score(a, b, float(score))

        f.close()


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by calling "align()"

    NOTE ON IMPLEMENTATION: You might find it helpful for code structure/organization to create additional functions
    that handle subtasks during the alignment. This is totally acceptable. However, you MUST implement the following
    functions with the expected behavior as outlined in the docstrings, which we will check in order to give
    at least partial credit for your P1 implementations:

    - populate_score_matrices
    - update_m, update_ix, and update_iy
    - find_traceback_start

    NOTE 2: Don't forget about that fuzzy_equals function at the top.

    """

    def __init__(self, input_file: str, output_file: str) -> None:
        """
        Initialize Align object.

        :param input_file: alignment input file path
        :param output_file: file path to write the output alignment
        """

        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters()

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        ### FILL IN ###
        # NOTE: be careful about how you initialize these!

        # Creating dimensions for
        self.nrow = len(self.align_params.seq_a) + 1
        self.ncol = len(self.align_params.seq_b) + 1

        self.m_matrix = ScoreMatrix('m', self.nrow, self.ncol)
        self.ix_matrix = ScoreMatrix('ix', self.nrow, self.ncol)
        self.iy_matrix = ScoreMatrix('iy', self.nrow, self.ncol)

    def align(self):
        """
        Main method for running the alignment.

        Note: there is no static typing on the method, as you can choose to return arbitrary
        intermediates/output if it's helpful for debugging. The essential minimal functionality that this
        method must have is to write the resulting alignments to the output file in the format specified in the
        project page
        """

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # print("-------M MATRIX-------")
        # self.m_matrix.print_scores()
        # # print("-------IX MATRIX-------")
        # self.ix_matrix.print_scores()
        # # print("-------IY MATRIX-------")
        # self.iy_matrix.print_scores()

        # print(self.find_traceback_start())

        # Perform traceback, then write the output to file
        self.traceback()
        self.write_output()

        # perform a traceback and write the output to an output file
        ### FILL IN ###

    def populate_score_matrices(self) -> None:
        """
        Populate the score matrices based on the data in align_params. Should call update(i,j) for each entry
        in the score matrices.
        """

        ### FILL IN ###

        # loop through from F(1,1) and call update() each time until populated
        for i in range(1, self.nrow):
            for j in range(1, self.ncol):
                self.update(i, j)

    def update(self, row: int, col: int) -> None:
        """
        Update all matrices at a given row and column index.

        :param row: index of row to update
        :param col: index of column to update
        """
        # Update each of the 3 matrices
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row: int, col: int) -> None:
        """
        Update matrix M.

        :param row: row index
        :param col: column index
        """

        # Extract the alignment score from the match matrix
        score = self.align_params.match_matrix.get_score(self.align_params.seq_a[row - 1],
                                                         self.align_params.seq_b[col - 1])

        # Store the values of the relevant prior squares in the m matrix
        m_prev = self.m_matrix.get_score(row - 1, col - 1)
        ix_prev = self.ix_matrix.get_score(row - 1, col - 1)
        iy_prev = self.iy_matrix.get_score(row - 1, col - 1)

        # Take the maximum value, according to our algorithm for the m matrix update function
        maximum = max(m_prev + score, ix_prev + score, iy_prev + score)

        # If local, the maximum should not be negative
        if not self.align_params.global_alignment and maximum < 0:
            maximum = 0

        # Update the m matrix with the maximum value
        self.m_matrix.set_score(row, col, maximum)

        # Record the pointers--where the maximum value came from
        pointers = set()
        if fuzzy_equals(maximum, m_prev + score):
            pointers.add((row - 1, col - 1, "m"))
        if fuzzy_equals(maximum, ix_prev + score):
            pointers.add((row - 1, col - 1, "ix"))
        if fuzzy_equals(maximum, iy_prev + score):
            pointers.add((row - 1, col - 1, "iy"))

        self.m_matrix.set_pointers(row, col, pointers)

    def update_ix(self, row: int, col: int) -> None:
        """
        Update matrix Ix.

        :param row: row index
        :param col: column index
        """

        # Store the values of the relevant prior squares in the m matrix
        m_prev = self.m_matrix.get_score(row - 1, col) - self.align_params.dy
        ix_prev = self.ix_matrix.get_score(row - 1, col) - self.align_params.ey

        # Take the maximum value, according to our algorithm for the ix matrix update function
        maximum = max(m_prev, ix_prev)

        # If local, the maximum should not be negative
        if not self.align_params.global_alignment and maximum < 0:
            maximum = 0

        # Update the matrix with the maximum value
        self.ix_matrix.set_score(row, col, maximum)

        # Record the pointers--where the maximum value came from
        pointers = set()
        if fuzzy_equals(maximum, m_prev):
            pointers.add((row - 1, col, "m"))
        if fuzzy_equals(maximum, ix_prev):
            pointers.add((row - 1, col, "ix"))
        self.ix_matrix.set_pointers(row, col, pointers)

    def update_iy(self, row: int, col: int) -> None:
        """
        Update matrix Iy.

        :param row: row index
        :param col: column index
        """

        # Store the values of the relevant prior squares in the m matrix
        m_prev = self.m_matrix.get_score(row, col - 1) - self.align_params.dx
        iy_prev = self.iy_matrix.get_score(row, col - 1) - self.align_params.ex

        # Take the maximum value, according to our algorithm for the ix matrix update function
        maximum = max(m_prev, iy_prev)

        # If local, the maximum should not be negative
        if not self.align_params.global_alignment and maximum < 0:
            maximum = 0

        # Update the matrix with the maximum value
        self.iy_matrix.set_score(row, col, maximum)

        # Record the pointers--where the maximum value came from
        pointers = set()
        if fuzzy_equals(maximum, m_prev):
            pointers.add((row, col - 1, "m"))
        if fuzzy_equals(maximum, iy_prev):
            pointers.add((row, col - 1, "iy"))
        self.iy_matrix.set_pointers(row, col, pointers)

    def find_traceback_start(self) -> Tuple[float, Set[Tuple[int, int, str]]]:
        """
        Find the location(s) to start the traceback and the corresponding best score.
        NOTE: Think carefully about how to set this up for local alignment.

        :return: The value of the best score and the location(s) to start the traceback to produce this score.
        The expected format is (max_val, max_loc), where max_val is the best score and max_loc is a set containing
        the positions to start the traceback that produce the best score, where each position is represented by a tuple
        [ex. (5.5, {(1,2), (3,4)}) ].
        """

        # If local, initiate the maximum to equal zero and the path to be (0, 0)
        if not self.align_params.global_alignment:
            maximum = 0
            tracks = {(0, 0, 'm')}

            # Locating the right traceback starting location. Loop through all the values in the matrix.
            for col in range(self.ncol):
                for row in range(self.nrow):
                    # If the value at a location is larger than the current match, update the max
                    if self.m_matrix.get_score(row, col) > maximum:
                        maximum = self.m_matrix.get_score(row, col)
                        # replace the traceback start with the new location
                        tracks = {(row, col, 'm')}
                    # If one of the values is equal to the max, simply add that location to the set
                    elif fuzzy_equals(self.m_matrix.get_score(row, col), maximum):
                        tracks.add((row, col, 'm'))
            print(maximum, tracks)
            return maximum, tracks

        # If Global, initiate the max to be the score at one of the valid locations in the m matrix
        else:
            maximum = self.m_matrix.get_score(self.nrow - 1, self.ncol - 1)
            tracks = {(self.nrow - 1, self.ncol - 1, 'm')}

            # Locating the right traceback starting location. Loop through all the values in the matrix.
            for col in range(self.ncol):
                for row in range(self.nrow):
                    # Only look through the last column and last row of the matrices
                    if col != (self.ncol - 1) and row != (self.nrow - 1):
                        continue
                    # Check the scores at valid locations in each of the matrices
                    for (matrix, name) in [(self.m_matrix, "m"), (self.ix_matrix, "ix"), (self.iy_matrix, "iy")]:
                        # If the value at a location is larger than the current match, update the max
                        if matrix.get_score(row, col) > maximum:
                            maximum = matrix.get_score(row, col)
                            tracks = {(row, col, name)}
                        # If one of the values is equal to the max, simply add that location to the set
                        elif fuzzy_equals(matrix.get_score(row, col), maximum):
                            tracks.add((row, col, name))
            print(maximum, tracks)
            return maximum, tracks

    def traceback_recurse(self, pointer: Tuple[int, int, str]):

        # Base case for Global: if we arrive at the first row or first column
        row, col, matrix = pointer
        if self.align_params.global_alignment and (row == 0 or col == 0):
            # Return an empty dict for each seq
            return [{'a': [], 'b': []}]

        # Update the scores based on the location of the pointer
        score = None
        if matrix == 'm':
            score = self.m_matrix.get_score(row, col)
        elif matrix == 'ix':
            score = self.ix_matrix.get_score(row, col)
        elif matrix == 'iy':
            score = self.iy_matrix.get_score(row, col)

        # Base case for Local: if the score is zero
        if not self.align_params.global_alignment and fuzzy_equals(score, 0):
            # Return an empty dict for each seq
            return [{'a': [], 'b': []}]

        # Initiating and updating the growing sequence:
        a = None
        b = None
        # If there is a match
        if matrix == 'm':
            a = self.align_params.seq_a[row - 1]
            b = self.align_params.seq_b[col - 1]
        # If there is a gap
        elif matrix == 'ix':
            a = self.align_params.seq_a[row - 1]
            b = '_'
        # If there is a gap
        elif matrix == 'iy':
            b = self.align_params.seq_b[col - 1]
            a = '_'

        # Updating the next pointer
        if matrix == "m":
            pointers = self.m_matrix.get_pointers(row, col)
        elif matrix == "ix":
            pointers = self.ix_matrix.get_pointers(row, col)
        elif matrix == "iy":
            pointers = self.iy_matrix.get_pointers(row, col)
        else:
            raise Exception("!")

        # Perform the recursive step for all of the pointers
        all_sub_alignments = []
        for pointer in pointers:
            sub_alignments = self.traceback_recurse(pointer)
            # add the values to traceback paths individually once they've split
            for sa in sub_alignments:
                all_sub_alignments.append(sa)

        # Add values to traceback paths once they've joined
        for sa in all_sub_alignments:
            sa['a'].append(a)
            sa['b'].append(b)

        return all_sub_alignments

    def max_gap(self):

        end = 0
        gap = 0
        start = 0
        real_start = 0
        counts = 0
        total_skips = 0

        for seq in self.sequences:
            for i in range(len(seq)):
                if seq[i] == '\n':
                    break
                if seq[i] == '_':
                    total_skips += 1
                    if counts == 0:
                        start = i
                    counts += 1
                elif seq[i] != '_' and counts > gap:
                    gap = counts
                    end = i-1 - total_skips
                    real_start = start - total_skips
                    start = 0
                    counts = 0
            total_skips = 0
            gap = 0
            counts = 0

        print("Max gap location on bacteria: ", real_start, " : ", end)
        print("gap: ", gap, " start: ", start, " counts: ", counts)

    def traceback(self):  ### FILL IN additional arguments ###
        """
        Perform a traceback.

        NOTE: There is no static typing on the method, as you have freedom to choose how you'd like to implement
        this method, including which arguments to include and what to return.

        HINT: It is extremely helpful for debugging to include a way to print the traceback path.
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """

        maximum, tracks = self.find_traceback_start()

        all_seqs = []
        for track in tracks:
            subseqs = self.traceback_recurse(track)
            for ss in subseqs:
                all_seqs.append(ss)

        # put results in correct format for grading
        seqs = [(''.join(d['a']) + '\n' + ''.join(d['b'])) for d in all_seqs]

        # store the results in attributes
        self.sequences = set(seqs)
        self.max_score = maximum

        self.max_gap()

    def write_output(self) -> None:
        """
        Write the output of an alignment to the output file.
        """
        with open(self.output_file, 'w') as output:
            output.write(str(round(self.max_score, 1)) + '\n')

            for seq in self.sequences:
                output.write('\n' + seq + '\n')

# DO NOT MODIFY THE CODE BELOW!
def main():
    """
    Run the align function from command line, passing the input and output paths as arguments.
    """

    # Check that arguments are in present in the command line as expected
    if (len(sys.argv) != 3):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__ == "__main__":
    main()