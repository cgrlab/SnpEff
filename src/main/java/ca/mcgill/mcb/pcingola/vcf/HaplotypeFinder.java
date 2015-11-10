package ca.mcgill.mcb.pcingola.vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ca.mcgill.mcb.pcingola.genotypes.GenotypeVector;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A 'haplotype' is a collection of variants that we analyze together
 * to infer a 'haplotype' (a.k.a. 'compound') annotation
 *
 * @author pcingola
 */
public class HaplotypeFinder {

	public static boolean debug = true;

	Map<String, GenotypeVector> genotypeVectors; // Sample genotypes
	Set<Haplotype> haplotypes; // Haplotypes found so far

	public HaplotypeFinder() {
		genotypeVectors = new HashMap<>();
		haplotypes = new HashSet<>();
	}

	/**
	 * Add new genotype vector and calculate all new haplotypes
	 * @return A list of all new found haplotypes
	 */
	public List<Haplotype> add(GenotypeVector gv) {
		// First find all new haplotypes that can be created using this new entry
		List<Haplotype> newhaplosGvs = haplotypesFromGenotypeVector(gv);

		// Find all haplotypes based on current haplotypes 
		Map<Haplotype, Haplotype> newHapsByOldHaps = haplotypesFromCurrentHaplotypes(gv);

		// Add calculated haplotypes from GenotypeVectors
		for (Haplotype hap : newhaplosGvs) {
			haplotypes.add(hap);
			if (debug) Gpr.debug("Adding haplotype:" + hap);
		}

		// Replace 'current' haplotypes by new ones in 'haplotypes'
		for (Haplotype hapOld : newHapsByOldHaps.keySet()) {
			Haplotype hapNew = newHapsByOldHaps.get(hapOld);
			haplotypes.remove(hapOld);
			haplotypes.add(hapNew);
			if (debug) Gpr.debug("Replacing haplotypes:" //
					+ "\n\tOld haplotype: " + hapOld //
					+ "\n\tNew haplotype: " + hapNew //
			);
		}

		// Add new entry
		addGenotypeVector(gv);

		// Return all found haplotypes
		if (debug) Gpr.debug(this + "\n\n");
		newhaplosGvs.addAll(newHapsByOldHaps.values());
		return newhaplosGvs;
	}

	void addGenotypeVector(GenotypeVector gv) {
		if (debug) Gpr.debug("Adding genotype: '" + toString(gv.getVariant()) + "'\t" + gv);
		genotypeVectors.put(toString(gv.getVariant()), gv);
	}

	protected boolean filter(Variant var) {
		return true;
	}

	GenotypeVector getGenotypeVector(Variant variant) {
		return genotypeVectors.get(toString(variant));
	}

	public Set<Haplotype> getHaplotypes() {
		return haplotypes;
	}

	/**
	 * Find new haplotypes using genotypes from 'genotypeVector'
	 */
	List<Haplotype> haplotypesFromGenotypeVector(GenotypeVector genotypeVector) {
		List<Haplotype> haplotypes = new LinkedList<>();

		// Compare to all current genotypeVectors
		for (GenotypeVector gv : genotypeVectors.values()) {
			if (genotypeVector.hasHaplotype(gv)) {
				Haplotype h = new Haplotype();
				h.add(gv.getVariant());
				h.add(genotypeVector.getVariant());
				haplotypes.add(h);
			}
		}

		return haplotypes;
	}

	/**
	 * Find new haplotypes using current haplotypes
	 */
	Map<Haplotype, Haplotype> haplotypesFromCurrentHaplotypes(GenotypeVector genotypeVector) {
		Map<Haplotype, Haplotype> newHapByOldHap = new HashMap<>();

		// Compare to all current haplotypes
		for (Haplotype hapOld : haplotypes) {
			if (debug) Gpr.debug("Analyzing haplotye: " + hapOld);

			// Create a list of genotypes form all variants in 'hap'
			List<GenotypeVector> gvs = new ArrayList<>();
			for (Variant var : hapOld) {
				GenotypeVector gv = getGenotypeVector(var);
				if (gv == null) throw new RuntimeException("Could not find GenotypeVector for variant " + toString(var));
				gvs.add(getGenotypeVector(var));
			}

			// Now check if there is a new haplotype 
			if (genotypeVector.hasHaplotype(gvs)) {
				// Create new haplotype with all the variants
				Haplotype hapNew = new Haplotype();
				hapNew.add(genotypeVector.getVariant());
				for (GenotypeVector gv : gvs)
					hapNew.add(gv.getVariant());

				newHapByOldHap.put(hapOld, hapNew);
			}
		}

		return newHapByOldHap;
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

		sb.append("Genotype vectors: " + genotypeVectors.size() + "\n");
		int i = 0;
		for (GenotypeVector g : genotypeVectors.values()) {
			sb.append(i + "\t" + g + "\n");
			i++;
		}

		if (!haplotypes.isEmpty()) {
			sb.append("Haplotypes: " + haplotypes.size() + "\n");
			i = 0;
			for (Haplotype h : haplotypes) {
				sb.append(i + "\t" + h + "\n");
				i++;
			}
		}

		return sb.toString();
	}

	String toString(Variant var) {
		return var.getChromosomeName() //
				+ ":" + var.getStart()//
				+ "_" + var.getReference() //
				+ "/" + var.getAlt() //
				;
	}

}
