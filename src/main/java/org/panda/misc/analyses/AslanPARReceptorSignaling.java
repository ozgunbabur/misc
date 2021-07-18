package org.panda.misc.analyses;

import org.panda.causalpath.run.ViewNeighborhoodInPriors;
import org.panda.misc.causalpath.CausalPathPenaltyBoxLister;
import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class AslanPARReceptorSignaling
{
	public static final String DATA_DIR = "/Users/ozgun/Documents/Data/Aslan/Thrombin-PAR/";
	public static final String BASE = "/Users/ozgun/Documents/Analyses/Aslan-Thrombin-PAR/";

	static Set<String> idMemory = new HashSet<>();


	public static void main(String[] args) throws IOException
	{
//		convertData();
//		explorePriorNeighborhood();
		writePenalyBoxLists();
	}

	public static void convertData() throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(BASE + "data.csv");
		writer.write("ID\tSymbols\tSites\tModification\tEffect\tResting-vs-PAR1\tResting-vs-PAR4\tResting-vs-Thrombin");

		String inFile = DATA_DIR + "ASLA-952_Quantitative_Summary.csv";
		String[] topHeader = FileUtil.readHeader(inFile);
		int par1Start = ArrayUtil.indexOf(topHeader, "Resting vs PAR1");
		int par4Start = ArrayUtil.indexOf(topHeader, "Resting vs PAR4");
		int thrombinStart = ArrayUtil.indexOf(topHeader, "Resting vs Thrombin");

		String[] header = FileUtil.lines(inFile).skip(5).findFirst().get().split("\t");

		int siteIndex = ArrayUtil.indexOf(header, "Site List");
		int symIndex = ArrayUtil.indexOf(header, "UniProt Gene Name");
		int fdrPar1Index = ArrayUtil.indexOf(header, par1Start, "FDR");
		int fdrPar4Index = ArrayUtil.indexOf(header, par4Start, "FDR");
		int fdrThrombinIndex = ArrayUtil.indexOf(header, thrombinStart, "FDR");
		int fcPar1Index = ArrayUtil.indexOf(header, par1Start, "FC");
		int fcPar4Index = ArrayUtil.indexOf(header, par4Start, "FC");
		int fcThrombinIndex = ArrayUtil.indexOf(header, thrombinStart, "FC");

		FileUtil.linesTabbed(inFile).skip(6).forEach(t ->
		{
			String sym = t[symIndex].split(" ")[0];

			List<String> siteList = Arrays.asList(t[siteIndex].split("; "));
			String id = getID(sym, siteList);

			double pPar1 = Double.valueOf(t[fdrPar1Index]);
			double pPar4 = Double.valueOf(t[fdrPar4Index]);
			double pThrombin = Double.valueOf(t[fdrThrombinIndex]);

			if (t[fcPar1Index].startsWith("-")) pPar1 *= -1;
			if (t[fcPar4Index].startsWith("-")) pPar4 *= -1;
			if (t[fcThrombinIndex].startsWith("-")) pThrombin *= -1;

			FileUtil.lnwrite(id + "\t" + sym + "\t" + CollectionUtil.merge(siteList, "|") + "\tP\t", writer);
			FileUtil.tab_write(pPar1 + "\t" + pPar4 + "\t" + pThrombin, writer);
		});

		writer.close();
	}

	static String getID(String sym, List<String> sites)
	{
		String pref = sym + "-" + CollectionUtil.merge(sites, "-");

		String id = pref + "-P";
		int repCounter = 0;

		while (idMemory.contains(id))
		{
			id = pref + "-" + (++repCounter) + "-P";
		}

		idMemory.add(id);
		return id;
	}

	static void explorePriorNeighborhood() throws IOException
	{
		String priorFile = "/Users/ozgun/Documents/Temp/causal-priors.txt";
		String formatFile = BASE + "Resting-vs-Thrombin/causative.format";
		String gene = "MTOR";
		String outDir = BASE + "prior-views/" + gene;
		ViewNeighborhoodInPriors.main(new String[]{"neighborhood", priorFile, formatFile, outDir, gene});
	}

	static void writePenalyBoxLists() throws IOException
	{
		CausalPathPenaltyBoxLister.writePenaltyBoxAsList(BASE + "Resting-vs-PAR1");
		CausalPathPenaltyBoxLister.writePenaltyBoxAsList(BASE + "Resting-vs-PAR4");
		CausalPathPenaltyBoxLister.writePenaltyBoxAsList(BASE + "Resting-vs-Thrombin");
	}
}
