package ca.mcgill.mcb.pcingola.snpEffect.testCases.unity;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.genotypes.Genotypes;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.Haplotype;
import ca.mcgill.mcb.pcingola.vcf.HaplotypeFinder;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import junit.framework.Assert;

/**
 * Test cases for haplotype detection
 */
public class TestCasesHaplotype {

	boolean debug = false;
	boolean verbose = false;

	public TestCasesHaplotype() {
	}

	/**
	 * Check that a set of haplotypes are found
	 */
	public void checkHaplotypes(String vcfFile, Set<String> expectedAnns) {
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		HaplotypeFinder hf = new HaplotypeFinder();

		for (VcfEntry ve : vcf) {
			if (ve.isMultiallelic()) {
				if (verbose) Gpr.debug(ve);
				for (Variant var : ve.variants()) {
					String alt = var.getGenotype();
					Genotypes gv = new Genotypes(ve, var.getGenotype());
					if (verbose) Gpr.debug("\tAlt: '" + alt + "'\tGenotype vector: " + gv);
					hf.add(gv);
				}
			} else {
				Genotypes gv = new Genotypes(ve);
				if (verbose) Gpr.debug(ve + "\n\tGenotype vector: " + gv);
				hf.add(gv);
			}
		}

		// Are there any haplotypes?
		Assert.assertTrue("No haplotype found", hf.hasHaplotypes());

		// How many haplotypes have been found?
		Set<Haplotype> haplotypes = hf.getHaplotypes();

		Assert.assertEquals("One haplotype expected:" //
				+ "\nExpected :\n" + setToStr(expectedAnns) //
				+ "\nFound    :\n" + setToStr(haplotypes) //
				, expectedAnns.size(), haplotypes.size());

		// Check that haplotypes match
		Assert.assertEquals("Haplotypes do not match expected", setToStr(expectedAnns), setToStr(haplotypes));
	}

	/**
	 * Check for a single haplotype
	 */
	public void checkHaplotypes(String vcfFile, String expectedAnnGt) {
		HashSet<String> expectedAnns = new HashSet<>();
		expectedAnns.add(expectedAnnGt);
		checkHaplotypes(vcfFile, expectedAnns);
	}

	/**
	 * Check that no haplotype is found
	 */
	public void checkNoHaplotypes(String vcfFile) {
		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		HaplotypeFinder hf = new HaplotypeFinder();

		for (VcfEntry ve : vcf) {
			Genotypes gv = new Genotypes(ve);
			if (verbose) Gpr.debug(ve + "\n\tGenotype vector: " + gv);
			hf.add(gv);
		}

		Assert.assertFalse("Haplotype found (none expected)", hf.hasHaplotypes());

		Set<Haplotype> haplotypes = hf.getHaplotypes();
		Assert.assertEquals("One haplotype expected", 0, haplotypes.size());
	}

	/**
	 * Convert a set to a list of sorted strings
	 */
	@SuppressWarnings("rawtypes")
	String setToStr(Set objects) {
		ArrayList<String> list = new ArrayList<>();

		for (Object o : objects)
			list.add(o.toString());
		Collections.sort(list);

		StringBuilder sb = new StringBuilder();
		for (String s : list)
			sb.append("\t" + s + "\n");

		return sb.toString();
	}

	/**
	 * Two consecutive variants affecting the same codon
	 * Implicit phasing: both variants are homozygous ALT in one sample.
	 */
	@Test
	public void test_01() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_01.vcf";
		checkHaplotypes(vcfFile, "7:116596741_T>A + 7:116596742_G>C");
	}

	/**
	 * Two consecutive variants affecting the same codon
	 * Implicit phasing: First variant is homozygous ALT and second variant is heterozygous ALT
	 */
	@Test
	public void test_02() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_02.vcf";
		checkHaplotypes(vcfFile, "7:116596741_T>A + 7:116596742_G>C");
	}

	/**
	 * Two consecutive variants affecting the same codon
	 * Explicit phasing
	 */
	@Test
	public void test_03() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_03.vcf";
		checkHaplotypes(vcfFile, "7:116596741_T>A + 7:116596742_G>C");
	}

	/**
	 * Two consecutive variants affecting the same codon
	 * Implicit phasing: Chromosome is haploid (i.e. genotype calls a "1" instead of "1/1")
	 */
	@Test
	public void test_04() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_04.vcf";
		checkHaplotypes(vcfFile, "7:116596741_T>A + 7:116596742_G>C");
	}

	/**
	 * Two consecutive variants affecting the same codon (SNP + MNP)
	 * Explicit phasing: Variants in different chromosomes (maternal and paternal)
	 */
	@Test
	public void test_05() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_05.vcf";
		checkNoHaplotypes(vcfFile);
	}

	/**
	 * Two consecutive variants affecting the same codon  (SNP + MNP)
	 * Implicit phasing: All variants are homozygous ALT
	 */
	@Test
	public void test_06() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_06.vcf";
		checkHaplotypes(vcfFile, "7:116596741_T>A + 7:116596742_GT>AA");
	}

	/**
	 * Three consecutive variants affecting the same codon  (SNP + SNP + SNP)
	 * Implicit phasing: All variants are homozygous ALT
	 */
	@Test
	public void test_07() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_07.vcf";
		checkHaplotypes(vcfFile, "7:116596741_T>A + 7:116596742_G>A + 7:116596743_T>A");
	}

	/**
	 * Two consecutive multiallelic variants affecting the same codon
	 * Implicit phasing: All variants are homozygous ALT
	 */
	@Test
	public void test_08() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_08.vcf";
		Set<String> expected = new HashSet<>();
		expected.add("7:116596741_T>A + 7:116596742_G>C");
		expected.add("7:116596741_T>C + 7:116596742_G>C");

		checkHaplotypes(vcfFile, expected);
	}

	/**
	 * Three consecutive variants affecting the same codon  (SNP + SNP + SNP)
	 * Implicit phasing: All variants are homozygous ALT
	 * Multiple samples with different genotypes create different haplotypes
	 */
	@Test
	public void test_09() {
		Gpr.debug("Test");
		String vcfFile = "tests/haplotype_09.vcf";

		Set<String> expected = new HashSet<>();
		expected.add("7:116596741_T>A + 7:116596742_G>A + 7:116596743_T>A");
		expected.add("7:116596741_T>A + 7:116596743_T>A");
		expected.add("7:116596742_G>A + 7:116596743_T>A");
		checkHaplotypes(vcfFile, expected);
	}

}
