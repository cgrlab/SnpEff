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
 * @author pcingola
 */
public class HaplotypeFinder {

	public static boolean debug = false;

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
		if (newGtByOldGt.isEmpty()) return null;

		// Add or replace genotypes
		boolean addGt = true;
		for (Genotypes gtOld : newGtByOldGt.keySet()) {
			Genotypes gtNew = newGtByOldGt.get(gtOld);

			// We can remove the old haplotype if the genotypes match and all variants are included
			if (gtNew.getHaplotype().containsAll(gtOld.getHaplotype()) //
					&& gtNew.equalsGenotypes(gtOld) //
			) {
				if (debug) Gpr.debug("Replacing:" //
						+ "\n\tOld: " + gtOld //
						+ "\n\tNew: " + gtNew //
				);
				remove(gtOld);
			}

			// Add genotypes
			haplotypes.add(gtNew);

			// If one of gtNew is equal, then we should not add 'gt' by itself
			// (because it will be redundant)
			addGt &= !gt.equalsGenotypes(gtNew);
		}

		// Add new entry
		if (addGt) genotypesByVariant.put(gt.getHaplotype().toString(), gt);

		// Return all found haplotypes
		if (debug) Gpr.debug(this + "\n\n");
		Set<Haplotype> haps = new HashSet<>();
		for (Genotypes g : newGtByOldGt.values())
			haps.add(g.getHaplotype());

		return haps;
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

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("Genotype vectors: " + genotypesByVariant.size() + "\n");
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
