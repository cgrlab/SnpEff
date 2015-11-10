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

	Map<String, Genotypes> genotypesByHaplotype; // All current genotypes
	Set<Genotypes> haplotypes; // All current genotypes having haplotypes

	public HaplotypeFinder() {
		genotypesByHaplotype = new HashMap<>();
		haplotypes = new HashSet<>();
	}

	/**
	 * Add new genotype vector and calculate all new haplotypes
	 * @return A list of all new found haplotypes
	 */
	public Collection<Haplotype> add(Genotypes gt) {
		// Find new genotypes based on current ones
		Map<Genotypes, Genotypes> newGtByOldGt = haplotypes(gt);

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
				// Remove form both collections
				haplotypes.remove(gtOld);
				genotypesByHaplotype.remove(gtOld.getHaplotype().toString());
			}

			haplotypes.add(gtNew);

			// If one of gtNew is equal, then we should not add 'gt' by itself 
			// (because it will be redundant)
			addGt &= !gt.equalsGenotypes(gtNew);
		}

		// Add new entry
		if (addGt) genotypesByHaplotype.put(gt.getHaplotype().toString(), gt);

		// Return all found haplotypes
		if (debug) Gpr.debug(this + "\n\n");
		Set<Haplotype> haps = new HashSet<>();
		for (Genotypes g : newGtByOldGt.values())
			haps.add(g.getHaplotype());
		return haps;
	}

	//	void addgenotypes(Genotypes gt) {
	//		if (debug) Gpr.debug("Adding genotype: '" + gt.getHaplotype() + "'\t" + gt);
	//		genotypesByHaplotype.put(gt.getHaplotype().toString(), gt);
	//	}

	protected boolean filter(Variant var) {
		return true;
	}

	//	Genotypes getgenotypes(Haplotype hap) {
	//		return genotypesByHaplotype.get(hap.toString());
	//	}

	public Set<Haplotype> getHaplotypes() {
		Set<Haplotype> gts = new HashSet<>();
		for (Genotypes gt : haplotypes)
			gts.add(gt.getHaplotype());
		return gts;
	}

	/**
	 * Find new haplotypes using current haplotypes
	 */
	Map<Genotypes, Genotypes> haplotypes(Genotypes genotypes) {
		Map<Genotypes, Genotypes> newGtByOldGt = new HashMap<>();

		// Create a list of all current genotypes
		List<Genotypes> gts = new LinkedList<>();
		gts.addAll(genotypesByHaplotype.values());
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

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("Genotype vectors: " + genotypesByHaplotype.size() + "\n");
		int i = 0;
		for (Genotypes gt : genotypesByHaplotype.values()) {
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
}
