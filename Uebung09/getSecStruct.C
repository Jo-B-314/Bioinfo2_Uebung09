#include "getSecStruct.h"
#include <vector>
#include <map>
#include <tuple>
#include <BALL/STRUCTURE/secondaryStructureProcessor.h>
#include <BALL/KERNEL/residue.h>
#include <BALL/SYSTEM/fileSystem.h>
#include <BALL/SYSTEM/directory.h>
#include <boost/filesystem.hpp>

using namespace std;
using namespace BALL;
using namespace boost::filesystem;

/**
 * our main can be called ./getSecStruct -p path
 * this argument is not optional
 * using -h or --help you'll get help
 **/
int main(int argc, char *argv[])
{
	//we didn't write any tests again, because it is only a main function and we didn't know how to test
	//and we didn't see any reason to, because we see the results and see that it works if we call it
	string h = "-h";
	string help = "--help";
	for (int i = 0; i < argc; i++)
	{
		if (h.compare(argv[i]) == 0 || help.compare(argv[i]) == 0)
		{
			//-h or --help found
			cout << "Help window:\n-p[PATH] : required, path to a directory with .pdb files\n";
			cout << "-h/--help: open this help window\n";
			return 0;
		}
	}
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Use -h or --help to open the help menu.\n";
		return 1;
	}
	string p = "-p";
	if (p.compare(argv[1]))
	{
		cerr << "No known flag. Use -h or --help to open the help menu.\n";
		return 1;
	}
	/*ofstream output;
	output.open("../data.txt");
	output << "xi-6\t"
		   << "xi-5\t"
		   << "xi-4\t"
		   << "xi-3\t"
		   << "xi-2\t"
		   << "xi-1\t"
		   << "xi\t";
	output << "xi+1\t"
		   << "xi+2\t"
		   << "xi+3\t"
		   << "xi+4\t"
		   << "xi+5\t"
		   << "xi+6\t"
		   << "si\t"
		   << "svm\n";*/
	//we need our last column to train our svm in the R program
	path pfad(argv[2]);
	directory_iterator end_dir_it;
	for (directory_iterator dir_it(pfad); dir_it != end_dir_it; dir_it++)
	{
		if (is_regular_file(dir_it->path()))
		{
			PDBFile input(dir_it->path().string());
			System system;
			input >> system;

			//we open the fragment data base
			FragmentDB fragment_db("");

			//normalization of atom names
			system.apply(fragment_db.normalize_names);
			system.apply(fragment_db.add_hydrogens);
			system.apply(fragment_db.build_bonds);

			map<size_t, tuple<char, char>> type_map;
			if (system.getProtein(0))
			{
				//we want to iterate over every amino acid
				Protein* protein = system.getProtein(0);
				auto res_it = protein->beginResidue();
				map<size_t, tuple<char, int>> type_map;
				size_t map_length = 0;
				for(; + res_it; res_it ++) {
					if (res_it->isAminoAcid()) {
						//find contacts of amino acids 
						//(bond between both Calpha < 7 A)
						//map one letter amino acid with index for matrix
						auto name = Peptides::OneLetterCode(res_it->getName());
						//count bonds
						int bonds = 0;
						std::tuple<char, int> tuple(name, bonds);
						type_map[map_length] = tuple;
					}
					map_length++;

				}
			}
		}
	}

	return 0;
}
