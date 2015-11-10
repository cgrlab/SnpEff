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

	List<Variant> variants;

	public Haplotype() {
		variants = new ArrayList<Variant>();
	}

	public void add(Haplotype haplotype) {
		variants.addAll(haplotype.getVariants());
	}

	public void add(Variant var) {
		variants.add(var);
	}

	/**
	 * Apply all variants in this haplotype to marker
	 */
	public Marker apply(Marker marker) {
		// Variants must be applied in reverse order (i.e. from
		// the end of the marker to the begining). This is because
		// an InDel can change the marker coordinates.
		Collections.sort(variants, Collections.reverseOrder());

		for (Variant var : variants)
			if (marker != null) { // A large deletion can delete the whole marker
				marker = marker.apply(var);
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
		Collections.sort(variants, Collections.reverseOrder());

		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for (Variant var : variants) {
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
		Collections.sort(variants);

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
