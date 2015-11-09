package ca.mcgill.mcb.pcingola.vcf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Variant;

/**
 * A 'haplotype' is a collection of variants that we analyze together
 * to infer a 'haplotype' (a.k.a. 'compound') annotation
 *
 * @author pcingola
 */
public class Haplotype {

	List<Variant> variants;
	Map<Variant, VcfEntry> vcfEntryByVariant;

	public Haplotype() {
		variants = new ArrayList<Variant>();
		vcfEntryByVariant = new HashMap<>();
	}

	public void add(Variant var, VcfEntry ve) {
		variants.add(var);
		vcfEntryByVariant.put(var, ve);
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

	public List<Variant> getHaplotype() {
		return variants;
	}

	/**
	 * Get VcfEntry associated with variant 'var'
	 */
	public VcfEntry getVcfEntry(Variant var) {
		return vcfEntryByVariant.get(var);
	}

}
