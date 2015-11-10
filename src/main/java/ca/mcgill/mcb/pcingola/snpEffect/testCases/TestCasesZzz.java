package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.List;

import org.junit.Test;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.genotypes.GenotypeVector;
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
	boolean verbose = false;

	public TestCasesZzz() {
	}

	public void checkHaplotypes(String vcfFile, String expectedAnnGt) {
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		HaplotypeFinder hf = new HaplotypeFinder();

		for (VcfEntry ve : vcf) {
			GenotypeVector gv = new GenotypeVector(ve);
			if (verbose) Gpr.debug(ve + "\n\tGenotype vector: " + gv);
			hf.add(gv);
		}

		Assert.assertTrue("No haplotype found", hf.hasHaplotype());

		List<Haplotype> haplotypes = hf.haplotypes();
		Assert.assertEquals("One haplotype expected", 1, haplotypes.size());

		Haplotype haplotype = haplotypes.get(0);
		Assert.assertEquals("Haplotype ", expectedAnnGt, haplotype.getAnnGenotype());
	}

	public void checkNoHaplotypes(String vcfFile, String expectedAnnGt) {
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

	/**
	 * Two consecutive multiallelic variants affecting the same codon 
	 * Implicit phasing: All variants are homozygous ALT
	 */
	@Test
	public void test_08() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_08.vcf";
		checkHaplotypes(vcfFile, "C-7:116596741_T>A");
		throw new RuntimeException("INCORRECT CHECKING");
	}

}
