package ca.mcgill.mcb.pcingola.vcf;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ca.mcgill.mcb.pcingola.genotypes.GenotypeVector;
import ca.mcgill.mcb.pcingola.interval.Variant;

/**
 * A 'haplotype' is a collection of variants that we analyze together
 * to infer a 'haplotype' (a.k.a. 'compound') annotation
 *
 * @author pcingola
 */
public class HaplotypeFinder {

	List<GenotypeVector> genotypeVectors;

	public HaplotypeFinder() {
		genotypeVectors = new ArrayList<>();
	}

	public void add(GenotypeVector gv) {
		genotypeVectors.add(gv);
	}

	protected boolean filter(Variant var) {
		return true;
	}

	public List<Haplotype> haplotypes() {
		List<Haplotype> haplotypes = new LinkedList<>();

		for (int i = 0; i < genotypeVectors.size(); i++) {
			GenotypeVector gvi = genotypeVectors.get(i);
			for (int j = i + 1; j < genotypeVectors.size(); j++) {
				GenotypeVector gvj = genotypeVectors.get(j);
				if (gvi.hasHaplotype(gvj)) {
					Haplotype h = new Haplotype();
					h.add(gvi.getVariant());
					h.add(gvj.getVariant());
					haplotypes.add(h);
				}
			}
		}

		return haplotypes;
	}

	/**
	 * Is there any haplotype ?
	 */
	public boolean hasHaplotype() {
		for (int i = 0; i < genotypeVectors.size(); i++) {
			GenotypeVector gvi = genotypeVectors.get(i);
			for (int j = i + 1; j < genotypeVectors.size(); j++) {
				GenotypeVector gvj = genotypeVectors.get(j);
				if (gvi.hasHaplotype(gvj)) return true;
			}
		}

		return false;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Genotype vectors: " + genotypeVectors.size() + "\n");
		for (int i = 0; i < genotypeVectors.size(); i++)
			sb.append(i + "\t" + genotypeVectors.get(i) + "\n");
		return sb.toString();
	}

}
