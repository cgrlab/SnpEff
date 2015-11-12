package ca.mcgill.mcb.pcingola;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Variant;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;

public class Zzz extends SnpEff {

	public static void main(String[] args) {
		String vcfFile = Gpr.HOME + "/t2d1/vcf/13k_clean.eff.vcf.gz";

		CountByType countType = new CountByType();
		CountByType countEff = new CountByType();
		int multiallelic = 0;
		int lineNum = 0;

		String chrPrev = "";
		VcfEntry vePrev = null;
		int endPrev = -1;

		VcfFileIterator vcf = new VcfFileIterator(vcfFile);
		for (VcfEntry ve : vcf) {

			if (ve.isMultiallelic()) multiallelic++;

			int end = 0;
			for (Variant var : ve.variants()) {
				countType.inc(var.getVariantType().toString());
				end = Math.max(end, var.getEnd());

				int diff = var.getStart() - endPrev;
				if (diff < 2) {

					int size = vcf.getSampleNames().size();
					int countHap = 0;
					for (int i = 0; i < size; i++) {
						VcfGenotype gt = ve.getVcfGenotype(i);
						VcfGenotype gtPrev = vePrev.getVcfGenotype(i);
						if ((gt.isHomozygousAlt() && gtPrev.isHeterozygous()) //
								|| (gtPrev.isHomozygousAlt() && gt.isHeterozygous())) {
							countHap++;
						}

					}

					// May be they are in the same codon?
					if (countHap > 0) Gpr.debug("SAME CODON: " + diff //
							+ "\tcount hap: " + countHap + " / " + size //
							+ "\n\t" + vePrev.toStringNoGt() //
							+ "\n\t" + ve.toStringNoGt() //
					);
				}
			}

			for (VcfEffect veff : ve.getVcfEffects()) {
				countEff.inc(veff.getEffectsStr());
			}

			if (lineNum > 1000) break;
			lineNum++;

			chrPrev = ve.getChromosomeName();
			vePrev = ve;
			endPrev = end;
		}

		System.out.println("Multi-allelic: " + multiallelic);
		System.out.println("Count by type:\n" + countType);
		System.out.println("Count by effect:\n" + countEff);
	}
}
