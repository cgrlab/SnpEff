package ca.mcgill.mcb.pcingola.vcf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Variant;

/**
 * A 'haplotype' is a collection of variants that we analyze together
 * to infer a 'haplotype' (a.k.a. 'compound') annotation
 *
 * @author pcingola
 */
public class Haplotype implements Iterable<Variant> {

	List<Variant> variants; // Variants sorted by start coordinate

	public Haplotype() {
		variants = new ArrayList<Variant>();
	}

	public void add(Haplotype haplotype) {
		variants.addAll(haplotype.getVariants());
	}

	public void add(Variant var) {
		variants.add(var);
		Collections.sort(variants); // Keep collection sorted
	}

	/**
	 * Apply all variants in this haplotype to marker
	 */
	public Marker apply(Marker marker) {
		// Variants must be applied in reverse order (i.e. from
		// the end of the marker to the begining). This is because
		// an InDel can change the marker coordinates.
		for (int i = size() - 1; i >= 0; i--) {
			if (marker != null) { // A large deletion can delete the whole marker
				Variant var = variants.get(i);
				marker = marker.apply(var);
			}
		}

		return marker;
	}

	/**
	 * Does this haplotype contain all variants in 'hap'?
	 */
	public boolean containsAll(Haplotype hap) {
		return variants.containsAll(hap.variants);
	}

	/**
	 * Genotype string in 'ANN' format
	 */
	public String getAnnGenotype() {
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for (int i = size() - 1; i >= 0; i--) {
			Variant var = variants.get(i);

			if (first) {
				sb.append(var.getGenotype());
				first = false;
			} else {
				sb.append("-" + var.getChromosomeName() //
						+ ":" + (var.getStart() + 1) //
						+ "_" + var.getReference() //
						+ ">" + var.getAlt()) //
						;
			}
		}

		return sb.toString();
	}

	/**
	 * Get start coordinate form "first" variant (i.e. the minimum start
	 * coordinate of any variant in this haplotype)
	 */
	public int getFirstStart() {
		int start = Integer.MAX_VALUE;
		for (Variant v : this)
			start = Math.min(start, v.getStart());
		return start;
	}

	/**
	 * Get start coordinate form "last" variant (i.e. the maximum start
	 * coordinate of any variant in this haplotype)
	 */
	public int getLastStart() {
		int start = 0;
		for (Variant v : this)
			start = Math.max(start, v.getStart());
		return start;
	}

	public List<Variant> getVariants() {
		return variants;
	}

	@Override
	public Iterator<Variant> iterator() {
		return variants.iterator();
	}

	public int size() {
		return variants.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Variant var : variants) {
			if (sb.length() > 0) sb.append(" + ");
			sb.append(var.getChromosomeName() //
					+ ":" + (var.getStart() + 1) //
					+ "_" + var.getReference() //
					+ ">" + var.getAlt()) //
					;
		}

		return sb.toString();
	}

}
