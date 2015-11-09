package ca.mcgill.mcb.pcingola.snpEffect.testCases;

import java.util.List;

import org.junit.Test;

import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEffCmdEff;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
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

	/**
	 * Annotate a VCF file and check that corresponding annotataion
	 */
	List<VcfEntry> annotate(String genome, String vcfFile) {
		String args[] = { "-noLog", "-noStats", genome, vcfFile };
		SnpEff snpEff = new SnpEff(args);
		snpEff.setVerbose(verbose);
		snpEff.setSupressOutput(!verbose);
		snpEff.setDebug(debug);

		SnpEffCmdEff seff = (SnpEffCmdEff) snpEff.snpEffCmd();
		List<VcfEntry> vcfEntries = seff.run(true);

		Assert.assertTrue("Empty annotataions list!", !vcfEntries.isEmpty());

		if (verbose) {
			Gpr.debug("Results:");
			for (VcfEntry ve : vcfEntries) {
				System.err.println("\t" + ve);
				for (VcfEffect eff : ve.getVcfEffects()) {
					System.err.println("\t\t" + eff.getHgvsP() + "\t" + eff);
				}
			}
		}

		return vcfEntries;
	}

	/**
	 * Two consecutive variants affecting the same codon
	 * Implicit phasing: both variants are homozygous ALT in one sample.
	 */
	@Test
	public void test_01() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_01.vcf";
		String genome = "test_ENST00000597499";
		List<VcfEntry> res = annotate(genome, vcfFile);

		// Check that the second entry is annotated as 'Cys => Thr'
		VcfEntry ve2 = res.get(1);
		VcfEffect veff2 = ve2.getVcfEffects().get(0);
		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Thr", veff2.getHgvsP());
	}

	//	/**
	//	 * Two consecutive variants affecting the same codon
	//	 * Implicit phasing: First variant is homozygous ALT and second variant is heterozygous ALT
	//	 */
	//	@Test
	//	public void test_02() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_02.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry is annotated as 'Cys => Thr'
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Thr", veff2.getHgvsP());
	//	}
	//
	//	/**
	//	 * Two consecutive variants affecting the same codon
	//	 * Explicit phasing
	//	 */
	//	@Test
	//	public void test_03() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_03.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry is annotated as 'Cys => Thr'
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Thr", veff2.getHgvsP());
	//	}
	//
	//	/**
	//	 * Two consecutive variants affecting the same codon
	//	 * Implicit phasing: Chromosome is haploid (i.e. genotype calls a "1" instead of "1/1")
	//	 */
	//	@Test
	//	public void test_04() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_04.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry is annotated as 'Cys => Thr'
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Thr", veff2.getHgvsP());
	//	}
	//
	//	/**
	//	 * Two consecutive variants affecting the same codon (SNP + MNP)
	//	 * Explicit phasing: Variants in different chromosomes
	//	 */
	//	@Test
	//	public void test_05() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_05.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry 
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Ser", veff2.getHgvsP());
	//	}
	//
	//	/**
	//	 * Two consecutive variants affecting the same codon  (SNP + MNP)
	//	 * Implicit phasing: All variants are homozygous ALT
	//	 */
	//	@Test
	//	public void test_06() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_06.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry 
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Lys", veff2.getHgvsP());
	//	}
	//
	//	/**
	//	 * Three consecutive variants affecting the same codon  (SNP + SNP + SNP)
	//	 * Implicit phasing: All variants are homozygous ALT
	//	 */
	//	@Test
	//	public void test_07() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_07.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry 
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Asn", veff2.getHgvsP());
	//
	//		// Check that the third entry 
	//		VcfEntry ve3 = res.get(1);
	//		VcfEffect veff3 = ve3.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Lys", veff2.getHgvsP());
	//	}

	//	/**
	//	 * Two consecutive multiallelic variants affecting the same codon 
	//	 * Implicit phasing: All variants are homozygous ALT
	//	 */
	//	@Test
	//	public void test_08() {
	//		Gpr.debug("Test");
	//		String vcfFile = "tests/haplotype_08.vcf";
	//		String genome = "test_ENST00000597499";
	//		List<VcfEntry> res = annotate(genome, vcfFile);
	//
	//		// Check that the second entry 
	//		VcfEntry ve2 = res.get(1);
	//		VcfEffect veff2 = ve2.getVcfEffects().get(0);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Thr", veff2.getHgvsP());
	//
	//		// Check that the second entry 
	//		VcfEffect veff2b = ve2.getVcfEffects().get(1);
	//		Assert.assertEquals("Incorrect haplotype annotation", "p.Cys13Pro", veff2b.getHgvsP());
	//	}

}
