package org.panda.misc;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.resource.tcga.MutTuple;
import org.panda.resource.tcga.MutationReader;
import org.panda.utility.statistics.BoxPlot;

import java.io.IOException;
import java.util.*;

/**
 * Generates the plot describing the effect of mutations to the expression.
 *
 * @author Ozgun Babur
 */
public class MutationToExpression
{
	public static final String TCGA_DIR = "/home/babur/Documents/TCGA/";
	public static final String STUDY = "BRCA";
	public static final String MUT_GENE = "TP53";
	public static final String EXP_GENE = "CDKN1A";
	public static final String OUT_DIR = "/home/babur/Documents/Temp/MutOnExp/";

	public static void main(String[] args) throws IOException
	{
		String dir = TCGA_DIR + STUDY + "/";
		MutationReader mr = new MutationReader(dir + "mutation.maf");
		Set<String> sampleSet = mr.getSamples();
		ExpressionReader er = new ExpressionReader(dir + "expression.txt");
		sampleSet.retainAll(er.getSamples());

		List<String> sampleList = new ArrayList<>(sampleSet);
		Collections.sort(sampleList);
		String[] samples = sampleList.toArray(new String[sampleList.size()]);
//		boolean[] muts = mr.getGeneAlterationArray(MUT_GENE, samples);
		boolean[] muts = new boolean[samples.length];

		List<MutTuple>[] mutations = mr.getMutations(MUT_GENE, samples);

		for (int i = 0; i < mutations.length; i++)
		{
			for (MutTuple tuple : mutations[i])
			{
				if (tuple.isDeleterious()) muts[i] = true;
			}
		}

		double[] exp = er.getGeneAlterationArray(EXP_GENE, samples);

		List<Double>[] valLists = new List[]{new ArrayList<>(), new ArrayList<>()};

		for (int i = 0; i < muts.length; i++)
		{
			if (muts[i]) valLists[0].add(exp[i]);
			else valLists[1].add(exp[i]);
		}

		BoxPlot.write(OUT_DIR + MUT_GENE + "-" + EXP_GENE + ".txt", new String[]{MUT_GENE + " mutated", MUT_GENE + " intact"}, valLists);
	}
}
