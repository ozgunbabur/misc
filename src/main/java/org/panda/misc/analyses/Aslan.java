package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Histogram;

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
	public static final String IN_BASE = "/home/ozgun/Data/Aslan-platelet/";

	public static final String OUT_BASE = "/home/ozgun/Documents/Analyses/Aslan-platelet/";

	public static final String[] INS_WITHOUT_FEEDBACK = new String[]{IN_BASE + "AslanTMT2017-GPVI-TiO2.csv", IN_BASE + "AslanTMT2017-GPVI-pTyr.csv"};
	public static final String[] INS_WITH_FEEDBACK = new String[]{IN_BASE + "AslanTMT2018-GPVIwithFeedback-TiO2.csv", IN_BASE + "AslanTMT2018-GPVIwithFeedback-pTyr.csv"};

	public static final String OUT_FILE_NAME = "/data-fdr0.1.txt";

	static final double FDR_THR = 0.1;

	public static void main(String[] args) throws IOException
	{
//		checkPhosphoProbDist();
		convertPhosphoproteomics();
//		printUniProtNames();
	}

	static void convertPhosphoproteomics() throws IOException
	{
		Map<String, Holder> map = new HashMap<>();
		for (String inF : INS_WITHOUT_FEEDBACK)
		{
			System.out.println("inF = " + inF);
			readInputWithoutFeedback(inF, map, FDR_THR);
		}

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(OUT_BASE + "without-feedback" + OUT_FILE_NAME));
		writer1.write("ID\tSymbols\tSites\tEffect\tFold change");

		map.values().stream().distinct().forEach(h -> FileUtil.lnwrite(h.toString(), writer1));

		writer1.close();

		map = new HashMap<>();
		for (String inF : INS_WITH_FEEDBACK)
		{
			System.out.println("inF = " + inF);
			readInputWithFeedback(inF, map, FDR_THR);
		}

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(OUT_BASE + "with-feedback" + OUT_FILE_NAME));
		writer2.write("ID\tSymbols\tSites\tEffect\tFold change");

		map.values().stream().distinct().forEach(h -> FileUtil.lnwrite(h.toString(), writer2));

		writer2.close();
	}

	static void readInputWithoutFeedback(String file, Map<String, Holder> map, double fdrThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(file)).skip(2).findFirst().get().split("\t");

		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "fold change");
		int logfcInd = ArrayUtil.indexOf(header, "logFC");
		int protInd = ArrayUtil.indexOf(header, "Protein Group Accessions");
		int sequenceInd = ArrayUtil.indexOf(header, "Sequence");
		int phosphoRSInd = ArrayUtil.indexOf(header, "phosphoRS Site Probabilities");
		int siteInd = ArrayUtil.indexOf(header, "Site List", "protein mod site");

		Files.lines(Paths.get(file)).skip(3).map(l -> l.split("\t"))
			.forEach(t ->
		{
			String prots = t[protInd];
			String sequence = t[sequenceInd].toUpperCase();
			String probs = t[phosphoRSInd];

			if (probs.startsWith("Too")) return;

			List<Integer> pepSites = getPepSites(probs);

			if (pepSites.isEmpty())
			{
				System.out.println("No sites for " + prots);
				return;
			}

			List<String> syms = new ArrayList<>();
			List<List<String>> sites = new ArrayList<>();

			String[] tt = prots.split(";");
			for (int i = 0; i < tt.length; i++)
			{
				String prot = tt[i];
				fillGenesAndSites(prot, sequence, pepSites, syms, sites);
			}

			if (syms.isEmpty())
			{
				if (!t[siteInd].isEmpty())
				{
					System.out.println("Cannot resolve genes and sites, but there is a site : " + t[siteInd]);
					fillGenesAndSites(tt[0], sequence, pepSites, syms, sites);
				}
				return;
			}

			double fdr = Double.parseDouble(t[fdrInd]);

			double fc = Double.parseDouble(t[fcInd]);
			if (t[logfcInd].startsWith("-")) fc = -fc;

			if (fdr > fdrThr) fc = 0;

			Holder h = new Holder(syms, sites, fc);

			// check for conflicting data
			if (map.containsKey(h.id) && map.get(h.id).fc * h.fc < 0)
			{
				System.err.println("Conflicting data for " + h.id);
			}

			if (!map.containsKey(h.id) || Math.abs(map.get(h.id).fc) < Math.abs(h.fc))
			{
				map.put(h.id, h);
			}
		});
	}

	static void readInputWithFeedback(String file, Map<String, Holder> map, double fdrThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(file)).skip(3).findFirst().get().split("\t");

		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "fold change");
		int logfcInd = ArrayUtil.indexOf(header, "logFC");
		int protInd = ArrayUtil.indexOf(header, "Full Accessions");
		int sequenceInd = ArrayUtil.indexOf(header, "New Sequence");
		int phosphoRSInd = ArrayUtil.indexOf(header, "New Site Prob Peptide");
		int siteInd = ArrayUtil.indexOf(header, "Site List", "protein mod site");

		Files.lines(Paths.get(file)).skip(4).map(l -> l.split("\t"))
			.forEach(t ->
		{
			String prots = t[protInd];
			String sequence = t[sequenceInd].toUpperCase().replaceAll("\\*", "");
			String probs = t[phosphoRSInd];

//			if (probs.startsWith("Too")) return;

			List<Integer> pepSites = getPepSites(probs);

			if (pepSites.isEmpty())
			{
				System.out.println("No sites for " + prots);
				return;
			}

			List<String> syms = new ArrayList<>();
			List<List<String>> sites = new ArrayList<>();

			String[] tt = prots.split(";");
			for (int i = 0; i < tt.length; i++)
			{
				String prot = tt[i].split("\\|")[1];
				fillGenesAndSites(prot, sequence, pepSites, syms, sites);
			}

			if (syms.isEmpty())
			{
				if (!t[siteInd].isEmpty())
				{
					System.out.println("Cannot resolve genes and sites, but there is a site : " + t[siteInd] + ", prot = " + prots);
				}
				return;
			}

			double fdr = Double.parseDouble(t[fdrInd]);

			double fc = Double.parseDouble(t[fcInd]);
			if (t[logfcInd].startsWith("-")) fc = -fc;

			if (fdr > fdrThr) fc = 0;

			Holder h = new Holder(syms, sites, fc);

			// check for conflicting data
			if (map.containsKey(h.id) && map.get(h.id).fc * h.fc < 0)
			{
				System.err.println("Conflicting data for " + h.id);
			}

			if (!map.containsKey(h.id) || Math.abs(map.get(h.id).fc) < Math.abs(h.fc))
			{
				map.put(h.id, h);
			}
		});
	}

	private static void fillGenesAndSites(String prot, String sequence, List<Integer> pepSites, List<String> syms, List<List<String>> sites)
	{
		String sym = HGNC.get().getSymbol(prot);
		if (sym != null)
		{

			int offset = UniProtSequence.get().getStartLocation(prot, sequence);

			if (offset < 0) System.out.println("Unknown UniProt ID = " + prot);

			if (offset > 0)
			{
				syms.add(sym);
				List<String> siteList = new ArrayList<>();
				for (Integer pepSite : pepSites)
				{
					String s = sequence.charAt(pepSite - 1) + "" + (offset + pepSite - 1);
					siteList.add(s);
				}
				sites.add(siteList);
			}
		}
	}

	static List<Integer> getPepSites(String phosphoRS)
	{
		List<Integer> list = new ArrayList<>();
		String[] t = phosphoRS.split("; ");
		for (int i = 0; i < t.length; i++)
		{
			double val = Double.valueOf(t[i].split(":")[1].trim());
			if (val >= 40)
			{
				String sites = t[i].substring(t[i].indexOf("(") + 1, t[i].indexOf(")"));

				for (String s : sites.split(","))
				{
					int num = Integer.valueOf(s);
					list.add(num);
				}
			}
		}
		return list;
	}

	static class Holder
	{
		List<String> syms;
		String id;
		double fc;
		List<List<String>> sites;

		Holder(List<String> syms, List<List<String>> sites, double fc)
		{
			this.syms = syms;
			this.fc = fc;
			this.sites = sites;

			String id = "";
			for (int i = 0; i < syms.size(); i++)
			{
				if (!id.isEmpty()) id += "-";
				id += syms.get(i);
				for (String s : sites.get(i))
				{
					id += "_" + s;
				}
			}

			this.id = id;
		}

		@Override
		public String toString()
		{
			return id + "\t" + CollectionUtil.merge(syms, " ") + "\t" + getSiteString() + "\t\t" + fc;
		}

		private String getSiteString()
		{
			String s = "";
			for (List<String> siteList : sites)
			{
				if (!s.isEmpty()) s += " ";
				s += CollectionUtil.merge(siteList, "|");
			}
			return s;
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

	static void checkPhosphoProbDist() throws IOException
	{
		Histogram h = new Histogram(10);
		h.setBorderAtZero(true);
		Files.lines(Paths.get("/home/ozgun/Data/Aslan-platelet/AslanTMT2017-GPVI-TiO2.csv")).skip(6)
			.map(l -> l.split("\t")).filter(t -> t.length > 4 && !t[4].startsWith("Too"))
			.map(t -> t[4].split("; ")).forEach(t ->
		{
//			System.out.println("t = " + Arrays.toString(t));
			for (int i = 0; i < t.length; i++)
			{
				h.count(Double.valueOf(t[i].split(" ")[1]));
			}
		});

		h.print();
	}
}
