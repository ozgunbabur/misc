package org.panda.misc.analyses;

import org.panda.causalpath.run.CausalityAnalysisSingleMethodInterface;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Aslan
{
//	public static final String BASE = "/home/babur/Documents/Analyses/Aslan/second-pass/";
//	public static final String PHOS_IN1 = BASE + "AslanPhosphoSerThrTMT-August2017.csv";
//	public static final String PHOS_IN2 = BASE + "AslanPhosphoTyrosineTMT-June2017.csv";
//	public static final String PHOS_OUT = BASE + "phosphoproteomics-fdr0.2.txt";

//	public static final String OLD_BASE = "/home/babur/Documents/Analyses/Aslan/first-pass/";
//	public static final String OLD_PHOS_IN = OLD_BASE + "phosphoproteomics.csv";
//	public static final String OLD_PHOS_OUT = OLD_BASE + "phosphoproteomics-fdr0.1.txt";

	public static final String BASE = "/home/babur/Documents/Analyses/Aslan/platelet-first-redo/";
	public static final String PHOS_IN1 = BASE + "ASLA-515_humanplatelet_TiO2phos_TO-SEND.csv";
	public static final String PHOS_IN2 = BASE + "ASLA-515_humanplatelet_pTyr_phos_TO-SEND.csv";
	public static final String PHOS_OUT = BASE + "data-fdr0.1.txt";


	static final double FDR_THR = 0.1;

	public static void main(String[] args) throws IOException
	{
		convertPhosphoproteomics(FDR_THR, PHOS_OUT, PHOS_IN1, PHOS_IN2);
//		convertPhosphoproteomics(FDR_THR, PHOS_OUT, PHOS_IN1, PHOS_IN2);
//		convertPhosphoproteomics(FDR_THR, OLD_PHOS_OUT, OLD_PHOS_IN);
//		printUniProtNames();
	}

	static void convertPhosphoproteomics(double fdrThr, String outFile, String... inFile) throws IOException
	{
		Map<String, Holder> map = new HashMap<>();
		for (String inF : inFile)
		{
			readCVSFile(inF, map, fdrThr);
		}


		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tFold change");

		map.values().stream().distinct().forEach(h -> FileUtil.lnwrite(h.toString(), writer));

		writer.close();
	}

	static void readCVSFile(String file, Map<String, Holder> map, double fdrThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(file)).skip(4).findFirst().get().split("\t");

		int siteInd = ArrayUtil.indexOf(header, "protein mod site", "Site List", "Mod Protein Location");
		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "fold change");
		int logfcInd = ArrayUtil.indexOf(header, "logFC");
		int geneInd = ArrayUtil.indexOf(header, "Protein Descriptions");
		int secondGeneInd = ArrayUtil.indexOf(header, "Gene List");

		Files.lines(Paths.get(file)).skip(5).map(l -> l.split("\t"))
			.forEach(t ->
		{
			String sym = t[geneInd];

			if (sym.contains("GN="))
			{
				int start = sym.indexOf("GN=") + 3;
				sym = sym.substring(start, sym.indexOf(" ", start));
			}
			else
			{
				sym = t[secondGeneInd];
				if (sym.contains(";")) sym = sym.substring(0, sym.indexOf(";"));
				if (sym.contains(" ")) sym = sym.substring(0, sym.indexOf(" "));

				if (sym.startsWith("Synonyms=") || sym.contains(":"))
				{
					System.err.println("No gene name = " + sym);
					return;
				}
			}

			// there must be a site in this data
			if (t[siteInd].isEmpty()) return;

			String[] sites = t[siteInd].split("; ");

			// Fix for an error in the old proteomics files
//			for (int i = 0; i < sites.length; i++)
//			{
//				sites[i] = oneDown(sites[i]);
//			}

			double fdr = Double.parseDouble(t[fdrInd]);

			double fc = Double.parseDouble(t[fcInd]);
			if (t[logfcInd].startsWith("-")) fc = -fc;

			if (fdr > fdrThr) fc = 0;

			String id = sym;
			for (String site : sites) id += "_" + site;

			Holder h = new Holder(id, sym, fc);

			Collections.addAll(h.sites, sites);

			// check for conflicting data
			if (map.containsKey(id) && map.get(id).fc * h.fc < 0)
			{
				System.err.println("Conflicting data for " + id);
			}

			if (!map.containsKey(id) || Math.abs(map.get(id).fc) < Math.abs(h.fc))
			{
				map.put(id, h);
			}
		});
	}

	static class Holder
	{
		String sym;
		String id;
		double fc;
		List<String> sites;

		Holder(String id, String sym, double fc)
		{
			this.sym = sym;
			this.id = id;
			this.fc = fc;
			sites = new ArrayList<>();
		}

		@Override
		public String toString()
		{
			return id + "\t" + sym + "\t" + ArrayUtil.getString("|", sites.toArray(new String[sites.size()])) +
				"\t\t" + fc;
		}
	}

	static String oneDown(String site)
	{
		String aa = site.substring(0, 1);
		int num = Integer.valueOf(site.substring(1));
		num--;
		return aa + num;
	}

	static void printUniProtNames() throws IOException
	{
		Files.lines(Paths.get("/home/babur/Documents/Analyses/Platelet-paper-2017/from-paper.csv")).skip(5)
			.map(l -> l.split("\t")[1].split(" ")[0]).forEach(System.out::println);
	}
}
