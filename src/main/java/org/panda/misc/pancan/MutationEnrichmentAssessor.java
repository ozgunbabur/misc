package org.panda.misc.pancan;

import org.panda.resource.tcga.MutTuple;
import org.panda.resource.tcga.MutationReader;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 *
 * @author Ozgun Babur
 */
public class MutationEnrichmentAssessor
{
	public static Map<String, Double>[] getEnrichmentPvalues(String maffile) throws IOException
	{
		MutationReader mr = new MutationReader(maffile);
		String[] samples = mr.getSamples().toArray(new String[0]);

		int mCnt = 0;
		int dCnt = 0;

		Map<String, Integer> mMap = new HashMap<>();
		Map<String, Integer> dMap = new HashMap<>();

		Set<String> genes = mr.getGenes();
		Progress p = new Progress(genes.size(), "Calculating enrichment p-values");
		for (String gene : genes)
		{
			p.tick();

			List<MutTuple>[] muts = mr.getMutations(gene, samples);

			List<String> types = Arrays.stream(muts).filter(Objects::nonNull).flatMap(Collection::stream)
				.map(m -> m.type).collect(Collectors.toList());

			int d = CollectionUtil.countValues(types, "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site");
			int m = CollectionUtil.countValues(types, "Missense_Mutation", "In_Frame_Del", "In_Frame_Ins");

			dMap.put(gene, d);
			mMap.put(gene, m);
			mCnt += m;
			dCnt += d;
		}

		Map<String, Double> dP = new HashMap<>();
		Map<String, Double> mP = new HashMap<>();

		int size = mCnt + dCnt;

		for (String gene : mMap.keySet())
		{
			int selected = dMap.get(gene) + mMap.get(gene);
			double pval_D = FishersExactTest.calcEnrichmentPval(size, dCnt, selected, dMap.get(gene));
			double pval_M = FishersExactTest.calcEnrichmentPval(size, mCnt, selected, mMap.get(gene));

			if (gene.equals("MTOR"))
			{
				System.out.print("");
			}

			if (pval_D + pval_M < 0.5)
			{
				int d = dMap.get(gene);
				int c = selected - d;
				int b = dCnt - d;
				int a = size - b - c + d;

				double pval = ChiSquare.testDependence(new long[][]{{a, b}, {c, d}});

				double globalRat = dCnt / (double) size;
				double currentRat = dMap.get(gene) / (double) selected;

				if (currentRat > globalRat)
				{
					pval_D = pval;
					pval_M = 1;
				}
				else
				{
					pval_D = 1;
					pval_M = pval;
				}
			}

			dP.put(gene, pval_D);
			mP.put(gene, pval_M);
		}

		return new Map[]{dP, mP};
	}

	public static void main(String[] args) throws IOException
	{
//		checkStatisticalUniformity();
		Map<String, Double>[] pvals = getEnrichmentPvalues("/home/babur/Documents/TCGA/PanCan/mutation.maf");

		Map<String, Double> dP = pvals[0];
		Map<String, Double> mP = pvals[1];

		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/babur/Documents/PanCan/mut-type-enrichment.txt"));
		writer.write("Gene\tInh Mut Enrich P-val\tMiss sense enrich P-val");
		dP.keySet().stream().sorted((g1, g2) ->
			new Double(Math.min(dP.get(g1), mP.get(g1))).compareTo(Math.min(dP.get(g2), mP.get(g2))))
			.forEach(g -> FileUtil.lnwrite(g + "\t" + dP.get(g) + "\t" + mP.get(g), writer));
		writer.close();

		System.out.println("FDR.select(dP, null, 0.05).size() = " + FDR.select(dP, null, 0.05).size());
		System.out.println("FDR.select(mP, null, 0.05).size() = " + FDR.select(mP, null, 0.05).size());
	}
}
