from align import *
import sys

def main():
    """
    Run the align function from command line, passing the input and output paths as arguments.
    """

    # Check that arguments are in present in the command line as expected
    if (len(sys.argv) != 4):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    output_file1 = sys.argv[1]
    output_file2 = sys.argv[2]
    output_file3 = sys.argv[3]
    # output_file4 = sys.argv[4]
    # output_file5 = sys.argv[5]

    # create an align object and run
    with open(output_file1) as f:
        lines = f.readlines()
        print('TURTLE: ', lines[0])

    with open(output_file2) as f:
        lines = f.readlines()
        print('CHIMP: ', lines[0])

    with open(output_file3) as f:
        lines = f.readlines()
        print('MOUSE: ', lines[0])

    # with open(output_file4) as f:
    #     lines = f.readlines()
    #     print('TRYPSIN: ', lines[0])
    #
    # with open(output_file5) as f:
    #     lines = f.readlines()
    #     print('UBX: ', lines[0])


if __name__ == "__main__":
    main()
