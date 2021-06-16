#include "getSecStruct.h"
#include "Matrix.h"
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
	path pfad(argv[2]);
	directory_iterator end_dir_it;
	//we have 20 row for each AA and one row for the Unknown ones
	//we start with up to 199 possible contacts, if we find more we will resize
	int poss_contacts = 200;
	Matrix matrix(21, 200);
	//N is the number of all residues (we need it for b))
	int N = 0;
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

			if (system.getProtein(0))
			{
				//we want to iterate over every amino acid
				Protein *protein = system.getProtein(0);
				auto res_it = protein->beginResidue();
				for (; + res_it; res_it++)
				{
					if (res_it->isAminoAcid())
					{
						N++;
						/*
                         * Alanine   A  0
                         * Glycine   G  1
                         * Valine    V  2
                         * Isoleucine I 3
                         * Leucine   L  4
                         * Methionin M  5
                         * Prolin    P  6
                         * Serine    S  7
                         * Threonine T  8
                         * Cysteine  C  9
                         * Phenylalanine F 10
                         * Tyrosine  Y  11
                         * Tryptophan W 12
                         * Lysine    K  13
                         * Arginine  R  14
                         * Histedine H  15
                         * Aspartate D  16
                         * Glutamate E  17
                         * Aspargine N  18
                         * Glutamine Q  19
                         */
						auto a_iti = res_it->beginAtom();
						Atom* C_alpha;
						int AA_index = 20;
						for (; + a_iti; a_iti++)
						{
							if (a_iti->getName() == "CA")
							{
								C_alpha = &*a_iti;
							}
						}
						//which is the index of our actual amino acid
						switch ((char)Peptides::OneLetterCode(res_it->getName()))
						{
						case 'A':
							AA_index = 0;
							break;
						case 'G':
							AA_index = 1;
							break;
						case 'V':
							AA_index = 2;
							break;
						case 'I':
							AA_index = 3;
							break;
						case 'L':
							AA_index = 4;
							break;
						case 'M':
							AA_index = 5;
							break;
						case 'P':
							AA_index = 6;
							break;
						case 'S':
							AA_index = 7;
							break;
						case 'T':
							AA_index = 8;
							break;
						case 'C':
							AA_index = 9;
							break;
						case 'F':
							AA_index = 10;
							break;
						case 'Y':
							AA_index = 11;
							break;
						case 'W':
							AA_index = 12;
							break;
						case 'K':
							AA_index = 13;
							break;
						case 'R':
							AA_index = 14;
							break;
						case 'H':
							AA_index = 15;
							break;
						case 'D':
							AA_index = 16;
							break;
						case 'E':
							AA_index = 17;
							break;
						case 'N':
							AA_index = 18;
							break;
						case 'Q':
							AA_index = 19;
							break;
						default:
							AA_index = 20;
							break;
						}
						//now we know the index for our actual aa and we have to count all our connections
						int contact_counter = 0;
						for (auto res_iterator = protein->beginResidue(); + res_iterator; res_iterator++)
						{
							if (!res_iterator->isAminoAcid())
							{
								continue;
							}
							else
							{
								Atom* otherCA;
								for (auto a_it = res_iterator->beginAtom(); + a_it; a_it++)
								{
									if (a_it->getName() == "CA")
									{
										otherCA = &*a_it;
									}
								}
								if (otherCA->getPosition().getDistance(C_alpha->getPosition()) < 10000)
								{
									//we have a contact
									contact_counter++;
								}
							}
						}
						//we find definitley a contact with our residue itself but we dont want to count it
						contact_counter--;
						//now we have to increase the value in our matrix but only if we have less contacts than our max value
						if (contact_counter < poss_contacts)
						{
							int val = matrix.getValue(AA_index, contact_counter);
							matrix.setValue(AA_index, contact_counter, val + 1);
						}
						else
						{
							matrix.resize(contact_counter + 1);
							poss_contacts = contact_counter + 1;
							int val = matrix.getValue(AA_index, contact_counter);
							matrix.setValue(AA_index, contact_counter, val + 1);
						}
					}
				}
			}
		}
	}
	cout << " ";
	for (int j = 0; j < poss_contacts; j++) {
		cout << j << " ";
	}
	cout << endl;
	for (int i = 0; i < 21; i++)
	{
		cout << getAAName(i) << " ";
		for (int j = 0; j < poss_contacts; j++)
		{
			cout << matrix.getValue(i, j) << " ";
		}
		cout << endl;
	}
	//now a) is finished lets start with b)
	//for each AA we compute how often it was seen
	Matrix matrix2(21, poss_contacts);
	for (int i = 0; i < 21; i++) {
		int Na = 0;
		//Na = times we found this aa
		for (int j = 0; j < poss_contacts; j++) {
			Na += matrix.getValue(i,j);
		}
		//Nk = times we found an aa with j contacts
		for (int j = 0; j < poss_contacts; j++) {
			int Nk = 0;
			for (int i2 = 0; i2 < 21; i2++){
				Nk += matrix.getValue(i2, j);
			}
			if ((Na*Nk) != 0) {
				matrix2.setValue(i, j, ((float)matrix.getValue(i, j) * N)/(Na*Nk));
			} else {
				matrix2.setValue(i, j, ((float)matrix.getValue(i, j) * N));
			}
		}
	}
	cout << " ";
	for (int j = 0; j < poss_contacts; j++) {
		cout << j << " ";
	}
	cout << endl;
	for (int i = 0; i < 21; i++)
	{
		cout << getAAName(i) << " ";
		for (int j = 0; j < poss_contacts; j++)
		{
			cout << matrix2.getValue(i, j) << " ";
		}
		cout << endl;
	}
	//now we come to c)
	
	return 0;
}
