package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.genotypes.GenotypeVector;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.Haplotype;
import ca.mcgill.mcb.pcingola.vcf.HaplotypeFinder;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import junit.framework.Assert;

/**
 * Test case
 *
 */
public class TestCasesZzz {

	boolean debug = false;
	boolean verbose = true;

	public TestCasesZzz() {
	}

	public void checkHaplotypes(String vcfFile, String expectedAnnGt) {
		HashSet<String> expectedAnns = new HashSet<>();
		expectedAnns.add(expectedAnnGt);
		checkHaplotypes(vcfFile, expectedAnns);
	}

	public void checkHaplotypes(String vcfFile, Set<String> expectedAnns) {
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		HaplotypeFinder hf = new HaplotypeFinder();

		for (VcfEntry ve : vcf) {
			if (ve.isMultiallelic()) {
				if (verbose) Gpr.debug(ve);
				for (Variant var : ve.variants()) {
					String alt = var.getGenotype();
					GenotypeVector gv = new GenotypeVector(ve, var.getGenotype());
					if (verbose) Gpr.debug("\tAlt: '" + alt + "'\tGenotype vector: " + gv);
					hf.add(gv);
				}
			} else {
				GenotypeVector gv = new GenotypeVector(ve);
				if (verbose) Gpr.debug(ve + "\n\tGenotype vector: " + gv);
				hf.add(gv);
			}
		}

		// Are there any haplotypes?
		Assert.assertTrue("No haplotype found", hf.hasHaplotype());

		// How many haplotypes have been found?
		List<Haplotype> haplotypes = hf.haplotypes();
		Assert.assertEquals("One haplotype expected", expectedAnns.size(), haplotypes.size());

		// Check that haplotypes match
		Set<String> haplotypesStr = new HashSet<>();
		for (Haplotype haplotype : haplotypes)
			haplotypesStr.add(haplotype.getAnnGenotype());

		Assert.assertEquals("Haplotypes do not match expected", setToStr(expectedAnns), setToStr(haplotypesStr));
	}

	String setToStr(Set<String> strs) {
		ArrayList<String> list = new ArrayList<>();
		list.addAll(strs);
		Collections.sort(list);

		StringBuilder sb = new StringBuilder();
		for (String s : list)
			sb.append("\t" + s + "\n");

		return sb.toString();
	}

	public void checkNoHaplotypes(String vcfFile) {
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		HaplotypeFinder hf = new HaplotypeFinder();

		for (VcfEntry ve : vcf) {
			GenotypeVector gv = new GenotypeVector(ve);
			if (verbose) Gpr.debug(ve + "\n\tGenotype vector: " + gv);
			hf.add(gv);
		}

		Assert.assertFalse("Haplotype found (none expected)", hf.hasHaplotype());

		List<Haplotype> haplotypes = hf.haplotypes();
		Assert.assertEquals("One haplotype expected", 0, haplotypes.size());
	}

	//	/**
	//	 * Three consecutive variants affecting the same codon  (SNP + SNP + SNP)
	//	 * Implicit phasing: All variants are homozygous ALT
	//	 */
	//	@Test
	//	public void test_07() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_07.vcf";
	//		checkHaplotypes(vcfFile, "C-7:116596741_T>A");
	//		throw new RuntimeException("INCORRECT CHECKING");
	//	}

}
