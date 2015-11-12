package ca.mcgill.mcb.pcingola.vcf;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ca.mcgill.mcb.pcingola.genotypes.Genotypes;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A 'haplotype' is a collection of variants that we analyze together
 * to infer a 'haplotype' (a.k.a. 'compound') annotation
 *
 *
 * Partial phasing and GATK's Read-Backed Phasing (RBP):
 * 		Reference: http://gatkforums.broadinstitute.org/discussion/45/purpose-and-operation-of-read-backed-phasing
 *
 * 		Phasing in the VCF format:
 * 		i) The "|" symbol is used for each sample to indicate that each of the alleles
 * 		of the genotype in question derive from the same haplotype as each of the
 * 		alleles of the genotype of the same sample in the previous NON-FILTERED variant
 * 		record. That is, rows without FILTER=PASS are essentially ignored in the read-backed
 * 		phasing (RBP) algorithm.
 *
 * 		ii) Note that the first heterozygous genotype record in a pair of haplotypes will
 * 		necessarily have a "/" - otherwise, they would be the continuation of the preceding
 * 		haplotypes.
 *
 * 		iii) A homozygous genotype is always "appended" to the preceding haplotype. For example, any
 * 		0/0 or 1/1 record is always converted into 0|0 and 1|1.
 *
 * 		iv) RBP attempts to phase a heterozygous genotype relative the preceding HETEROZYGOUS
 * 		genotype for that sample. If there is sufficient read information to deduce the two
 * 		haplotypes (for that sample), then the current genotype is declared phased ("/" changed
 * 		to "|") and assigned a PQ that is proportional to the estimated Phred-scaled error rate.
 * 		All homozygous genotypes for that sample that lie in between the two heterozygous genotypes
 * 		are also assigned the same PQ value (and remain phased).
 *
 * 		v) If RBP cannot phase the heterozygous genotype, then the genotype remains with a "/", and no
 * 		PQ score is assigned. This site essentially starts a new section of haplotype for this sample.
 *
 * @author pcingola
 */
public class HaplotypeFinder {

	public static boolean debug = true;

	Map<String, Genotypes> genotypesByVariant; // All current genotypes from only one variant (not haplotypes)
	Set<Genotypes> haplotypes; // All current genotypes having haplotypes (consisting of more than one variant)

	public HaplotypeFinder() {
		genotypesByVariant = new HashMap<>();
		haplotypes = new HashSet<>();
	}

	/**
	 * Add new genotype vector and calculate all new haplotypes
	 * @return A list of all new found haplotypes
	 */
	public Collection<Haplotype> add(Genotypes gt) {
		update(gt);
		if (!filter(gt.getHaplotype())) return null;

		// Find new genotypes based on current ones
		Map<Genotypes, Genotypes> newGtByOldGt = haplotypes(gt);
		if (newGtByOldGt.isEmpty()) {
			addGenotype(gt);
			return null;
		}

		// Add or replace genotypes
		boolean addGt = true;
		for (Genotypes gtOld : newGtByOldGt.keySet()) {
			Genotypes gtNew = newGtByOldGt.get(gtOld);

			// We can replace the old haplotype if the genotypes match and all variants are included
			if (gtNew.getHaplotype().containsAll(gtOld.getHaplotype()) && gtNew.equalsGenotypes(gtOld)) {
				replace(gtOld, gtNew);
			} else addHaplotypes(gtNew); // Don't replace, just add

			// If one of gtNew is equal, then we should not add 'gt' by itself
			// (because it will be redundant)
			addGt &= !gt.equalsGenotypes(gtNew);
		}

		// Add new entry
		if (addGt) addGenotype(gt);

		// Return all found haplotypes
		if (debug) Gpr.debug(this + "\n\n");
		Set<Haplotype> haps = new HashSet<>();
		for (Genotypes g : newGtByOldGt.values())
			haps.add(g.getHaplotype());

		return haps;
	}

	void addGenotype(Genotypes gt) {
		genotypesByVariant.put(gt.getHaplotype().toString(), gt);
	}

	// Add genotypes to 'haplotypes' collection
	void addHaplotypes(Genotypes gtNew) {
		haplotypes.add(gtNew);
	}

	/**
	 * Filter all variants in haplotype
	 * @return true if all variants pass filter
	 */
	public boolean filter(Haplotype hap) {
		for (Variant var : hap)
			if (!filter(var)) return false;
		return true;
	}

	/**
	 * Apply a filter to a variant
	 * @return True is variant passes filter
	 */
	public boolean filter(Variant var) {
		return true;
	}

	/**
	 * Filter VCF entry
	 */
	public boolean filter(VcfEntry ve) {
		return ve.isFilterPass();
	}

	public Set<Haplotype> getHaplotypes() {
		Set<Haplotype> gts = new HashSet<>();
		for (Genotypes gt : haplotypes)
			gts.add(gt.getHaplotype());
		return gts;
	}

	/**
	 * Find new haplotypes using current haplotypes
	 */
	protected Map<Genotypes, Genotypes> haplotypes(Genotypes genotypes) {
		Map<Genotypes, Genotypes> newGtByOldGt = new HashMap<>();

		// Create a list of all current genotypes
		List<Genotypes> gts = new LinkedList<>();
		gts.addAll(genotypesByVariant.values());
		gts.addAll(haplotypes);

		// Compare to all current haplotypes
		for (Genotypes gtOld : gts) {
			if (debug) Gpr.debug("Analyzing genotype: " + gtOld);

			// Now check if there is a new haplotype
			Genotypes gtNew = genotypes.haplotype(gtOld);
			if (gtNew != null) newGtByOldGt.put(gtOld, gtNew);
		}

		return newGtByOldGt;
	}

	/**
	 * Is there any haplotype ?
	 */
	public boolean hasHaplotypes() {
		return !haplotypes.isEmpty();
	}

	/**
	 * Remove genotypes form all collections
	 */
	protected void remove(Genotypes gt) {
		haplotypes.remove(gt);
		genotypesByVariant.remove(gt.getHaplotype().toString());
	}

	void replace(Genotypes gtOld, Genotypes gtNew) {
		if (debug) Gpr.debug("Replacing:" //
				+ "\n\tOld: " + gtOld //
				+ "\n\tNew: " + gtNew //
		);
		remove(gtOld);
		addHaplotypes(gtNew); // Add genotypes
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("Genotypes: " + genotypesByVariant.size() + "\n");
		int i = 0;
		for (Genotypes gt : genotypesByVariant.values()) {
			sb.append(i + "\t" + gt + "\n");
			i++;
		}

		if (!haplotypes.isEmpty()) {
			sb.append("Haplotypes: " + haplotypes.size() + "\n");
			i = 0;
			for (Genotypes gt : haplotypes) {
				sb.append(i + "\t" + gt + "\n");
				i++;
			}
		}

		return sb.toString();
	}

	/**
	 * Update caches and remove old entries
	 */
	protected void update(Genotypes gt) {
	}
}
