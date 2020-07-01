import argparse as ap
import parmed

parser = ap.ArgumentParser(description='Convert Amber to Gromacs')
parser.add_argument('infile', type=str, help='amber filename (w/out the extension)')
parser.add_argument('-itp', action = 'store_true', help='whether to break out an itp file, or to write one stand-alone .top')
args = parser.parse_args()  
prefix = args.infile

print("--- Reading in file `{}` ---".format(prefix))

amber = parmed.load_file(prefix+'.parm7', xyz=prefix+'.rst7')
if args.itp:
    print("Saving a separate .itp file with the atom,bond types and parameters")
    amber.save(prefix+'.top', parameters=prefix+'.itp', overwrite=True)
else:
    amber.save(prefix+'.top', overwrite=True)
amber.save(prefix+'.gro', overwrite=True)

