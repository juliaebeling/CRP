#include "gtest/gtest.h"
#include "mpc/Candidate.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/PhotoDisintegration.h"

#include <fstream>

namespace mpc {

TEST(ElectronPairProduction, EnergyDecreasing) {
	// test if energy loss occurs for protons with energies from 1e15 - 1e23 eV
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1,1)); // proton
	std::vector<Candidate *> secondaries;

	ElectronPairProduction epp1(ElectronPairProduction::CMB);
	for (int i = 0; i<80; i++) {
		double E = pow(10, 15+i*0.1) * eV;
		candidate.current.setEnergy(E);
		epp1.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() <= E);
	}

	ElectronPairProduction epp2(ElectronPairProduction::IR);
	for (int i = 0; i<80; i++) {
		double E = pow(10, 15+i*0.1) * eV;
		candidate.current.setEnergy(E);
		epp2.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}

	ElectronPairProduction epp3(ElectronPairProduction::CMBIR);
	for (int i = 0; i<80; i++) {
		double E = pow(10, 15+i*0.1) * eV;
		candidate.current.setEnergy(E);
		epp3.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}
}

TEST(ElectronPairProduction, BelowEnergyTreshold) {
	// test if nothing happens below 1e15 eV
	ElectronPairProduction epp(ElectronPairProduction::CMB);

	Candidate candidate;
	candidate.current.setId(getNucleusId(1,1)); // proton
	std::vector<Candidate *> secondaries;

	double E = 1e14 * eV;
	candidate.current.setEnergy(E);
	epp.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.current.getEnergy(), E);
}

TEST(ElectronPairProduction, EnergyLossValues) {
	// test if energy loss corresponds to table
	std::vector<double> x;
	std::vector<double> y;
	std::ifstream infile("data/ElectronPairProduction/cmbir.txt");
	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				x.push_back(a * eV);
				y.push_back(b * eV / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();

	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1,1)); // proton
	std::vector<Candidate *> secondaries;

	ElectronPairProduction epp;
	for (int i = 0; i<x.size(); i++) {
		candidate.current.setEnergy(x[i]);
		epp.process(&candidate, secondaries);
		double dE = x[i] - candidate.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE, dE_table, 1e-12);
	}
}


TEST(NuclearDecay, Neutron) {
	// test beta- decay n -> p
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1,0));
	candidate.current.setEnergy(1 * EeV);
	std::vector<Candidate *> secondaries;
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(1,1));
}

TEST(NuclearDecay, Scandium44) {
	// test beta+ decay 44Sc -> 44Ca
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(44,21));
	candidate.current.setEnergy(1 * EeV);
	std::vector<Candidate *> secondaries;
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(44,20));
}

TEST(NuclearDecay, Chlorium30) {
	// test proton emission 30Cl -> 29S
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(30,17));
	candidate.current.setEnergy(1 * EeV);
	std::vector<Candidate *> secondaries;
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(29,16));
}

TEST(NuclearDecay, Neon33) {
	// test neutron emission 33Ne -> 32Ne
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(33,10));
	candidate.current.setEnergy(1 * EeV);
	std::vector<Candidate *> secondaries;
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(32,10));
}

TEST(NuclearDecay, LimitNextStep) {
	// test if Decay module limits the step size for neutrons
	Candidate candidate;
	candidate.setNextStep(10 * Mpc);
	candidate.current.setId(getNucleusId(1,0));
	candidate.current.setEnergy(10 * EeV);
	std::vector<Candidate *> secondaries;
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_TRUE(candidate.getNextStep() < 10 * Mpc);
}

TEST(PhotoDisintegration, Basic) {
	PhotoDisintegration m;
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc
