package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.utility.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class NolanPaper
{
	public static final String BASE = "/Users/ozgun/Documents/Analyses/NolanPaper/";


	public static void main(String[] args) throws IOException
	{
//		convertPhosphoproteomics(BASE + "phospho-data.csv", BASE + "data.txt");
		takeSubset();
	}

	private static void convertPhosphoproteomics(String inFile, String outFile) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tSignedP");

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");

		int symIndex = ArrayUtil.indexOf(header, "GeneNames");
		int siteIndex = ArrayUtil.indexOf(header, "Phosphorylation position (>90% localization)");
		int mergedFCIndex = ArrayUtil.indexOf(header, "iTRAQ and TMT Combined Fold-Change");
		int mergedPIndex = ArrayUtil.indexOf(header, "iTRAQ and TMT Fisher's Combined P-Value");
		int p1Index = ArrayUtil.indexOf(header, "iTRAQ Adjusted P-Value");
		int fc1Index = ArrayUtil.indexOf(header, "iTRAQ Log Fold-Change");
		int p2Index = ArrayUtil.indexOf(header, "TMT Adjusted P-Value");
		int fc2Index = ArrayUtil.indexOf(header, "TMT Log Fold-Change");

		TermCounter tc = new TermCounter();

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[symIndex];
			if (sym.contains(" ")) return;

			String up = HGNC.get().getUniProt(sym);
			if (up == null) return;

			List<String> siteList = new ArrayList<>();
			String[] sites = t[siteIndex].split(";");
			for (String siteStr : sites)
			{
				if (siteStr.startsWith(up + "_"))
				{
					String loc = siteStr.split("_")[1];
					if (!loc.equals("NANA")) siteList.add(loc);
				}
			}

			tc.addTerm("" + siteList.size());
			if (siteList.isEmpty()) return;

			int pIndex = !t[mergedPIndex].equals("na") ? mergedPIndex : !t[p1Index].equals("na") ? p1Index : p2Index;
			int fcIndex = pIndex == mergedPIndex ? mergedFCIndex : pIndex == p1Index ? fc1Index : fc2Index;


			String signedP = t[pIndex];
			if (Double.valueOf(signedP) == 0D) signedP = "0.0000000001";
			if (t[fcIndex].startsWith("-")) signedP = "-" + signedP;


			String id = sym;
			for (String s : siteList)
			{
				id += "_" + s;
			}

			String siteStr = CollectionUtil.merge(siteList, "|");

			FileUtil.lnwrite(id + "\t" + sym + "\t" + siteStr + "\t\t" + signedP, writer);
		});

		writer.close();
		tc.print();
	}

	private static void takeSubset() throws IOException
	{
		CausalPathSubnetwork.writeGOINeighForCompBased(BASE + "fdr0.1", new HashSet<>(Arrays.asList("AKT1", "IRS1")), StreamDirection.BOTHSTREAM, "AKT1-IRS1");
	}

	private static void temp()
	{
		String up = HGNC.get().getUniProt("PYGB");
		System.out.println("up = " + up);
	}
}
