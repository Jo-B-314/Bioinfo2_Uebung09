#ifndef GET_SEQUENCE_H
#define GET_SEQUENCE_H

#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/chain.h>
#include <BALL/KERNEL/protein.h>
#include <BALL/FORMAT/PDBFile.h>
#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/STRUCTURE/peptides.h>
#include <BALL/common.h>
#include <BALL/DATATYPE/options.h>
#include <BALL/CONCEPT/embeddable.h>
#include <BALL/KERNEL/selector.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/PTE.h>

#include <iostream>
#include <string>

using namespace std;
using namespace BALL;

int c (Protein* protein) {
	cout << "Number of chains in the protein: " << protein->countChains() << endl;
	int str = protein->countChains();
	return str;
}

void d (Protein* protein) {
	int counter = 1;
		for (ChainIterator ch_it = protein->beginChain(); +ch_it; ch_it++) {
			cout << "Number of atoms in chain " << counter << ": " << ch_it->countAtoms() << endl;
			counter++;
		}
}

void e (System system) {
	Protein* protein = system.getProtein(0);
		for (ChainIterator ch_it = protein->beginChain(); +ch_it; ch_it++) {
			Peptides::OneLetterAASequence aa;
			//we want the AA Sequence for the whole protein
			auto re_it = ch_it->beginResidue();
			for (; +re_it; re_it++) {
				if (re_it->isAminoAcid()) {
					cout << Peptides::OneLetterCode(re_it->getName()) << " ";
					aa.push_back(Peptides::OneLetterCode(re_it->getName()));
				}
			}
			//we only want to print the sequence in one letter code
		}
}

void fg (System* system) {
	//we want a System pointer, because we want to delete the molecules also for the system in main function
	Selector wat("residue(HOH) OR residue(WAT)");
		system->apply(wat);
		//count
		int atom_num = wat.getNumberOfSelectedAtoms();
		//H2O is composed by 3 atoms
		int mol_num = atom_num / 3;
		cout << "Number of water molecules: " << mol_num << endl;
		//extension g: remove water molecules
		system->removeSelected();
		cout << "water molecules removed" << endl;
		system->deselect();
}

void task_h (System system) {
	Protein* protein = system.getProtein(0);
	int counter = 1;
	for (auto ch_it = protein->beginChain(); +ch_it; ch_it++) {
		bool aa_present = false;
		auto re_it = ch_it->beginResidue();
			for (; +re_it; re_it++) {
				if (re_it->isAminoAcid()) {
					aa_present = true;
					break;
				}
			}
			//there is no amino acid in the chain
			if (! aa_present) {
				cout << "No amino acids in chain " << counter << endl;
				//we want to know the number of C, N, O and H atoms
				AtomIterator ait = ch_it->beginAtom();
				int hcounter = 0;
				int ocounter = 0;
				int ccounter = 0;
				int ncounter = 0;
				for (; +ait; ait++) {
					String sym = ait->getElement().getSymbol();
					String H ("H");
					String C ("C");
					String N ("N");
					String O ("O");
					if (sym == C) {
						ccounter++;
						continue;
					}
					if (sym == N) {
						ncounter++;
						continue;
					}
					if (sym == O) {
						ocounter++;
						continue;
					}
					if (sym == H) {
						hcounter++;
					}
				}
				cout << "Number of C atoms: " << ccounter << endl;
				cout << "Number of O atoms: " << ocounter << endl;
				cout << "Number of N atoms: " << ncounter << endl;
				cout << "Number of H atoms: " << hcounter << endl;
			}
		counter++;
	}
}

char getAAName(int i) {
	switch(i) {
		case 0: return 'A';
		case 1: return 'G';
		case 2: return 'V';
		case 3: return 'I';
		case 4: return 'L';
		case 5: return 'M';
		case 6: return 'P';
		case 7: return 'S';
		case 8: return 'T';
		case 9: return 'C';
		case 10: return 'F';
		case 11: return 'Y';
		case 12: return 'W';
		case 13: return 'K';
		case 14: return 'R';
		case 15: return 'H';
		case 16: return 'D';
		case 17: return 'E';
		case 18: return 'N';
		case 19: return 'Q';
		default: return 'U';
	}
}

#endif