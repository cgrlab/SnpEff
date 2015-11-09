package ca.mcgill.mcb.pcingola.snpEffect.testCases.integration;

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
public class TestCasesHaplotypeAnn {

	boolean debug = false;
	boolean verbose = true;

	public TestCasesHaplotypeAnn() {
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
	 * Implicit phasing: both variants are homozygous REF in one sample.
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

}
