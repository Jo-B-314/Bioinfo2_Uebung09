#include <gtest/gtest.h>

#include <BALL/KERNEL/protein.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/KERNEL/system.h>
#include <BALL/STRUCTURE/peptideBuilder.h>
#include <vector>

#include "getSecStruct.h"

using namespace BALL;
using namespace std;

TEST(getSecStruct, c)
{
	vector<Peptides::AminoAcidDescriptor> descriptor_sequence;
	Peptides::AminoAcidDescriptor* aad = new Peptides::AminoAcidDescriptor;
	aad->setAminoAcidType("GLY");
	descriptor_sequence.push_back(*aad);
	Peptides::PeptideBuilder* pb = new Peptides::PeptideBuilder(descriptor_sequence);
    FragmentDB fdb("");
    pb->setFragmentDB(&fdb);
	Protein* protein = pb->construct();

	int str = c(protein);
	ASSERT_EQ(1, str);
}

/*TEST(getSecStruct, e)
{
	vector<Peptides::AminoAcidDescriptor> descriptor_sequence;
	Peptides::AminoAcidDescriptor* aad = new Peptides::AminoAcidDescriptor;
	aad->setAminoAcidType("GLY");
	descriptor_sequence.push_back(*aad);
	Peptides::PeptideBuilder* pb = new Peptides::PeptideBuilder(descriptor_sequence);
    FragmentDB fdb("");
    pb->setFragmentDB(&fdb);
	Protein* protein = pb->construct();

	System system;
	system.add(protein);
	system.apply(fdb.normalize_names);
	system.apply(fdb.add_hydrogens);
	system.apply(fdb.build_bonds);

	String res = e(system);


	ASSERT_EQ(0, str.copmare("G", res));

}*/