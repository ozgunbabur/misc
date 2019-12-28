package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.resource.signednetwork.SignedType;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.PhosphoGraph;
import org.panda.utility.statistics.Histogram;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class Aslan
{
	public static final String IN_BASE = "/home/ozgun/Data/Aslan/platelet-phosphoproteomics/";

	public static final String OUT_BASE = "/home/ozgun/Analyses/Aslan-platelet/";

	public static final String[] INS_WITHOUT_FEEDBACK = new String[]{IN_BASE + "AslanTMT2017-GPVI-TiO2.csv", IN_BASE + "AslanTMT2017-GPVI-pTyr.csv"};
	public static final String[] INS_WITH_FEEDBACK = new String[]{IN_BASE + "AslanTMT2018-GPVIwithFeedback-TiO2.csv", IN_BASE + "AslanTMT2018-GPVIwithFeedback-pTyr.csv"};

	public static final String OUT_FILE_NAME = "/data-fdr0.1.txt";

	static final double FDR_THR = 0.1;

	public static void main(String[] args) throws IOException
	{
//		checkPhosphoProbDist();
//		convertPhosphoproteomics();
//		printUniProtNames();
//		convertYoungAndBold();
//		convertIndividualOlds();
//		convertIndividualYoungs();

//		printOtherPenaltied();
//		listDarkPhosphoproteome();

//		printStatisticsFromDataFile();
//		printStatisticsFromResultsFiles();

//		readValidationResults();
		compareResultsToValidations();
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

	static void convertYoungAndBold() throws IOException
	{
		String inDir = IN_BASE + "boldandyoung/";
		String outBase = OUT_BASE + "youngandbold/";

		String[] cases = new String[]{"MaleVsFemale", "OldFemaleVsYoungFemale", "OldMaleVsOldFemale",
			"OldMaleVsYoungMale", "OldVsYoung", "YoungMaleVsYoungFemale"};

		for (String aCase : cases)
		{
			Map<String, Holder> map = new HashMap<>();
			readPhosphoproteomicsOldAndYoung(inDir + aCase + "_TiO2phos_TMT.csv", map, 0.1);
			readTotalProteomicsOldAndYoung(inDir + aCase + "_totalProteome_TMT.csv", map, 0.1);

			FileUtil.mkdirs(outBase + aCase);
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(outBase + aCase + "/data.txt"));

			writer.write("ID\tSymbols\tSites\tEffect\tFold change");
			map.values().stream().distinct().forEach(h -> FileUtil.lnwrite(h.toString(), writer));
			writer.close();

			ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});

			BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(outBase + aCase + "/graph.format"));
			writer2.write("node\tall-nodes\tcolor\t255 255 255");
			writer2.write("\nnode\tall-nodes\tbordercolor\t0 0 0");
			map.values().stream().distinct().filter(h -> h.fc != 0).forEach(h ->
				FileUtil.write(h.getFormatString(vtc), writer2));
			writer2.close();
		}
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

	static void readPhosphoproteomicsOldAndYoung(String file, Map<String, Holder> map, double fdrThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(file)).skip(4).findFirst().get().split("\t");

		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "Fold Change");
		int logfcInd = ArrayUtil.indexOf(header, "logFC");
		int geneInd = ArrayUtil.indexOf(header, "protein_description");
		int siteInd = ArrayUtil.indexOf(header, "absolute_mods");

		Files.lines(Paths.get(file)).skip(5).map(l -> l.split("\t"))
			.forEach(t ->
		{
			String gene = t[geneInd];
			int gInd = gene.indexOf(" GN=");
			if (gInd < 0) return;
			gene = gene.substring(gInd + 4, gene.indexOf(" ", gInd + 4));

			List<String> syms = new ArrayList<>();
			List<List<String>> sites = new ArrayList<>();

			List<String> list = new ArrayList<>();
			String siteStr = t[siteInd];
			for (String s : siteStr.split("; "))
			{
				if (s.endsWith("(Phospho)")) list.add(s.substring(0, s.indexOf("(")));
			}

			syms.add(gene);
			sites.add(list);

			if (syms.isEmpty())
			{
				if (!t[siteInd].isEmpty())
				{
					System.out.println("Cannot resolve genes and sites, but there is a site : " + t[siteInd] + ", prot = " + gene);
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

	static void readTotalProteomicsOldAndYoung(String file, Map<String, Holder> map, double fdrThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(file)).skip(4).findFirst().get().split("\t");

		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fcInd = ArrayUtil.indexOf(header, "Fold Change");
		int logfcInd = ArrayUtil.indexOf(header, "logFC");
		int geneInd = ArrayUtil.indexOf(header, "Description");

		Files.lines(Paths.get(file)).skip(5).map(l -> l.split("\t")).filter(t -> t.length > 2).forEach(t ->
		{
			String gene = t[geneInd];
			int gInd = gene.indexOf(" GN=");
			if (gInd < 0) return;
			gene = gene.substring(gInd + 4, gene.indexOf(" ", gInd + 4));

			List<String> syms = new ArrayList<>();
			List<List<String>> sites = new ArrayList<>();

			syms.add(gene);

			if (syms.isEmpty())
			{
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
				if (!sites.isEmpty())
				{
					for (String s : sites.get(i))
					{
						id += "_" + s;
					}
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

		private String getFormatString(ValToColor vtc)
		{
			StringBuilder sb = new StringBuilder();

			for (int i = 0; i < syms.size(); i++)
			{
				if (sites.isEmpty())
				{
					sb.append("\nnode\t").append(syms.get(i)).append("\tcolor\t").append(vtc.getColorInString(fc));
				}
				else
				{
					SiteEffectCollective sec = new SiteEffectCollective();
					Integer effect = null;
					for (String site : sites.get(i))
					{
						effect = sec.getEffect(syms.get(i), site);
						if (effect != null && effect != 0) break;
					}
					sb.append("\nnode\t").append(syms.get(i)).append("\trppasite\t").append(id).append("|p|").append(vtc.getColorInString(fc)).append("|").append(effect != null && effect == 1 ? "0 180 20" : effect != null && effect == -1 ? "180 0 20" : "0 0 0").append("|").append(fc);
				}
			}

			return sb.toString();
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



	static void convertIndividualOlds() throws IOException
	{
		String inDir = "/home/ozgun/Data/Aslan/platelet-phosphoproteomics/boldandyoung/OldIndividuals/";
		String outDir = "/home/ozgun/Analyses/Aslan-platelet/youngandbold/OldIndividuals/";

		convertIndividual(inDir + "OM7.csv", outDir + "OM7/data.txt", 0.1);
		convertIndividual(inDir + "OM8.csv", outDir + "OM8/data.txt", 0.1);
		convertIndividual(inDir + "OM9.csv", outDir + "OM9/data.txt", 0.1);
		convertIndividual(inDir + "OF10.csv", outDir + "OF10/data.txt", 0.1);
		convertIndividual(inDir + "OF13.csv", outDir + "OF13/data.txt", 0.1);
	}

	static void convertIndividualYoungs() throws IOException
	{
		String inDir = "/home/ozgun/Data/Aslan/platelet-phosphoproteomics/boldandyoung/YoungIndividuals/";
		String outDir = "/home/ozgun/Analyses/Aslan-platelet/youngandbold/YoungIndividuals/";

		convertIndividual(inDir + "YM1.csv", outDir + "YM1/data.txt", 0.1);
		convertIndividual(inDir + "YM2.csv", outDir + "YM2/data.txt", 0.1);
		convertIndividual(inDir + "YM3.csv", outDir + "YM3/data.txt", 0.1);
		convertIndividual(inDir + "YF4.csv", outDir + "YF4/data.txt", 0.1);
		convertIndividual(inDir + "YF14.csv", outDir + "YF14/data.txt", 0.1);
	}

	static void convertIndividual(String inFile, String outFile, double fdrThr) throws IOException
	{
		System.out.println("inFile = " + inFile);
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tFold change");

		String[] header = Files.lines(Paths.get(inFile)).skip(4).findFirst().get().split("\t");
		int symInd = ArrayUtil.indexOf(header, "protein_description");
		int siteInd = ArrayUtil.indexOf(header, "absolute_mods");
		int foldInd = ArrayUtil.indexOf(header, "fold change", "fold change ");
		int signInd = ArrayUtil.indexOf(header, "Log2(A/B)");
		int fdrInd = ArrayUtil.indexOf(header, "BH_adjusted");

		Files.lines(Paths.get(inFile)).skip(5).map(l -> l.split("\t")).forEach(t ->
		{
			String symbol = t[symInd];
			int ind = symbol.indexOf(" GN=");
			if (ind < 0) return;
			symbol = symbol.substring(ind + 4, symbol.indexOf(" ", ind + 5));

			String[] s = t[siteInd].split("; ");
			List<String> sites = new ArrayList<>();
			for (String ss : s)
			{
				if (ss.contains("(Phospho)"))
				{
					ss = ss.replace("(Phospho)", "");
					sites.add(ss);
				}
			}

			double value = 0;

			double fdr = Double.valueOf(t[fdrInd]);
			if (fdr <= fdrThr)
			{
				value = Double.valueOf(t[foldInd]);
				if (!t[signInd].startsWith("-")) value *= -1;
			}

			String siteStr = CollectionUtil.merge(sites, "|");
			FileUtil.lnwrite(symbol + "_" + siteStr + "\t" + symbol + "\t" + siteStr + "\t\t" + value, writer);
		});
		writer.close();
	}

	/**
	 * To get a list of proteins significantly changing but not even mapping to causal priors.
	 */
	public static void printOtherPenaltied() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Aslan-platelet/without-feedback/";
		String resultFile = dir + "causative.sif";
		String dataFile = dir + "data-fdr0.1.txt";

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10},
			new Color[]{new Color(40, 80, 255), Color.WHITE, new Color(255, 80, 40)});

		Set<String> existing = new HashSet<>();

		Files.lines(Paths.get(resultFile)).map(l -> l.split("\t")).forEach(t ->
		{
			existing.add(t[0]);
			if (t.length > 2) existing.add(t[2]);
		});

		SiteEffectCollective sec = new SiteEffectCollective();

		Set<String> genes = new HashSet<>();

		Files.lines(Paths.get(dataFile)).skip(1).map(l -> l.split("\t"))
			.filter(t -> !t[4].equals("0.0") && CollectionUtil.intersectionEmpty(existing,
				new HashSet<>(Arrays.asList(t[1].split(" "))))).forEach(t ->
		{
			String[] g = t[1].split(" ");
			String[] s = t[2].split(" ");

			for (int i = 0; i < g.length; i++)
			{
				String gene = g[i];
				genes.add(gene);
				String sites = s[i];

				Integer effect = getSiteEffect(sec, gene, sites);

				String id = t[0];
				double v = Double.valueOf(t[4]);

				String borderC = effect == null || effect == 0 ? "50 50 50" : effect == 1 ? "0 180 20" : "180 0 20";

				System.out.println("node\t"+ gene + "\trppasite\t" + id + "|p|" + vtc.getColorInString(v) + "|" + borderC + "|" + v);
			}

		});

		System.out.println("\n\n\n");
		for (String gene : genes)
		{
			System.out.println(gene);
		}
	}

	/**
	 * To get a list of proteins significantly changing but not even mapping to causal priors.
	 */
	public static void listDarkPhosphoproteome() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Aslan-platelet/with-feedback/";
		String resultFile = dir + "causative.sif";
		String dataFile = dir + "data-fdr0.1.txt";

		Set<String> inGraph = new HashSet<>();

		Files.lines(Paths.get(resultFile)).map(l -> l.split("\t")).forEach(t ->
		{
			if (t.length > 2)
			{
				inGraph.add(t[0]);
				inGraph.add(t[2]);
			}
		});

		SiteEffectCollective sec = new SiteEffectCollective();

		Files.lines(Paths.get(dataFile)).skip(1).map(l -> l.split("\t"))
			.filter(t -> !t[4].equals("0.0") && CollectionUtil.intersectionEmpty(inGraph,
				new HashSet<>(Arrays.asList(t[1].split(" "))))).sorted(Comparator.comparing(o -> o[1])).forEach(t ->
		{
			String[] g = t[1].split(" ");
			String[] s = t[2].split(" ");

			Integer effect = null;
			for (int i = 0; i < g.length; i++)
			{
				String gene = g[i];
				String sites = s[i];

				if (effect == null || effect == 0) effect = getSiteEffect(sec, gene, sites);
				else break;
			}

			System.out.println(t[1] + "\t" + t[2] + "\t" + t[4] + "\t" + ((effect == null || effect == 0) ? "" : effect == 1 ? "a" : "i"));

		});
	}

	private static Integer getSiteEffect(SiteEffectCollective sec, String gene, String sites)
	{
		Integer effect = null;
		for (String site : sites.split("\\|"))
		{
			Integer e = sec.getEffect(gene, site);
			if (e != null && e != 0)
			{
				effect = e;
				break;
			}
		}

		return effect;
	}

	private static void printStatisticsFromDataFile() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Aslan-platelet/";
		int[] cnt = new int[4];
		int rows = 0;
		int upreg = 1;
		int downreg = 2;
		int overlapInd = 3;
		Set<String> sitesC1 = new HashSet<>();
		Set<String> sitesC2 = new HashSet<>();
		Set<String> upregSitesC1 = new HashSet<>();
		Set<String> upregSitesC2 = new HashSet<>();
		Set<String> downregSitesC1 = new HashSet<>();
		Set<String> downregSitesC2 = new HashSet<>();
		Set<String> protsC1 = new HashSet<>();
		Set<String> protsC2 = new HashSet<>();

//		System.out.println("without-feedback SIF:");
//		Set<String> woOld = readSIF(dir + "without-feedback-relax1aa/causative.sif");
//		Set<String> woNew = readSIF(dir + "without-feedback-relax1aa-copy/causative.sif");
//		CollectionUtil.printVennCounts(woOld, woNew);
//		woOld.removeAll(woNew);
//		woOld.forEach(System.out::println);
//		System.out.println();

		Files.lines(Paths.get(dir + "without-feedback-relax1aa/data-fdr0.1.txt")).skip(1)
			.map(l -> l.split("\t")).forEach(t ->
		{
			cnt[rows]++;
			int sign = 0;
			if (!t[4].equals("0.0"))
			{
				if (t[4].startsWith("-"))
				{
					cnt[downreg]++;
					sign = -1;
				}
				else
				{
					cnt[upreg]++;
					sign = 1;
				}
			}

			sitesC1.addAll(distributeSites(t[0]));
			if (sign == 1) upregSitesC1.addAll(distributeSites(t[0]));
			else if (sign == -1) downregSitesC1.addAll(distributeSites(t[0]));

			protsC1.addAll(Arrays.asList(t[1].split(" ")));
		});

		System.out.println("without-feedback row stats:");
		System.out.println("cnt[rows] = " + cnt[rows]);
		System.out.println("cnt[upreg] = " + cnt[upreg]);
		System.out.println("cnt[downreg] = " + cnt[downreg]);

		for (int i = 0; i < 3; i++)
		{
			cnt[i] = 0;
		}

//		System.out.println("\nwith-feedback SIF:");
//		Set<String> wOld = readSIF(dir + "with-feedback-relax1aa/causative.sif");
//		Set<String> wNew = readSIF(dir + "with-feedback-relax1aa-copy/causative.sif");
//		CollectionUtil.printVennCounts(wOld, wNew);
//		wOld.removeAll(wNew);
//		wOld.forEach(System.out::println);


		Files.lines(Paths.get(dir + "with-feedback-relax1aa/data-fdr0.1.txt")).skip(1)
			.map(l -> l.split("\t")).forEach(t ->
		{
			cnt[rows]++;
			int sign = 0;
			if (!t[4].equals("0.0"))
			{
				if (t[4].startsWith("-"))
				{
					cnt[downreg]++;
					sign = -1;
				}
				else
				{
					cnt[upreg]++;
					sign = 1;
				}
			}

			if (distributeSites(t[0]).stream().anyMatch(sitesC1::contains)) cnt[overlapInd]++;
			sitesC2.addAll(distributeSites(t[0]));
			if (sign == 1) upregSitesC2.addAll(distributeSites(t[0]));
			else if (sign == -1) downregSitesC2.addAll(distributeSites(t[0]));
			protsC2.addAll(Arrays.asList(t[1].split(" ")));
		});

		System.out.println("with-feedback row stats:");
		System.out.println("\ncnt[rows] = " + cnt[rows]);
		System.out.println("cnt[upreg] = " + cnt[upreg]);
		System.out.println("cnt[downreg] = " + cnt[downreg]);

		System.out.println("\nwithout-feedback sites = " + sitesC1.size());
		System.out.println("with-feedback sites = " + sitesC2.size());
		System.out.println();
		CollectionUtil.printVennCounts(sitesC1, sitesC2);

		System.out.println("\nwithout-feedback upreg sites = " + upregSitesC1.size());
		System.out.println("with-feedback upreg sites = " + upregSitesC2.size());
		System.out.println();
		CollectionUtil.printVennCounts(upregSitesC1, upregSitesC2);

		System.out.println("\nwithout-feedback downreg sites = " + downregSitesC1.size());
		System.out.println("with-feedback downreg sites = " + downregSitesC2.size());
		System.out.println();
		CollectionUtil.printVennCounts(downregSitesC1, downregSitesC2);

		System.out.println("\noverlapping rows = " + cnt[overlapInd]);

		System.out.println("\nwithout-feedback prots = " + protsC1.size());
		System.out.println("with-feedback prots = " + protsC2.size());
		CollectionUtil.printVennCounts(protsC1, protsC2);
	}


	private static void printStatisticsFromResultsFiles() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Aslan-platelet/";

		Set<String> relsC1 = readSIF(dir + "without-feedback-relax1aa/causative.sif");
		Set<String> relsC2 = readSIF(dir + "with-feedback-relax1aa/causative.sif");
		System.out.println("relsC1.size() = " + relsC1.size());
		System.out.println("relsC2.size() = " + relsC2.size());
		CollectionUtil.printVennCounts(relsC1, relsC2);

		String file = dir + "without-feedback-relax1aa-copy/results.txt";
		Set<String> protsC1 = new HashSet<>();
		Set<String> sitesC1 = new HashSet<>();
		readResultFile(file, protsC1, sitesC1);

		file = dir + "with-feedback-relax1aa-copy/results.txt";
		Set<String> protsC2 = new HashSet<>();
		Set<String> sitesC2 = new HashSet<>();
		readResultFile(file, protsC2, sitesC2);

		System.out.println("\nprotsC1.size() = " + protsC1.size());
		System.out.println("protsC2.size() = " + protsC2.size());
		CollectionUtil.printVennCounts(protsC1, protsC2);
		
		System.out.println("\nsitesC1.size() = " + sitesC1.size());
		System.out.println("sitesC2.size() = " + sitesC2.size());
		CollectionUtil.printVennCounts(sitesC1, sitesC2);
		
	}

	private static void readResultFile(String file, Set<String> prots, Set<String> sites) throws IOException
	{
		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			prots.add(t[0]);
			prots.add(t[2]);

//			Set<String> relSites = relaxSites(t[3]);

			String sID = selectRelatedID(t[4], t[0]);
			String tID = selectRelatedID(t[7], t[2]);
			
			sites.addAll(distributeSites(sID));
			sites.addAll(distributeSites(tID));
		});
	}

	private static Set<String> readSIF(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).map(l -> l.split("\t")).filter(t -> t.length > 2).map(t -> ArrayUtil.getString("\t", t[0], t[1], t[2])).collect(Collectors.toSet());
	}

	private static Set<String> distributeSites(String id)
	{
		Set<String> set = new HashSet<>();
		String[] t = id.split("_");
		for (int i = 1; i < t.length; i++)
		{
			set.add(t[0] + "_" + t[i]);
		}
		return set;
	}

	private static Set<String> relaxSites(String relStr)
	{
		Set<String> relaxed = new HashSet<>();
		for (String site : relStr.split(";"))
		{
			String aa = site.substring(0, 1);
			int pos = Integer.parseInt(site.substring(1));
			relaxed.add(site);
			relaxed.add(aa + (pos - 1));
			relaxed.add(aa + (pos + 1));
		}
		return relaxed;
	}
	
	private static String selectRelatedID(String id, String sym)
	{
		if (id.endsWith("-by-network-sig") || id.endsWith("-active")) return id;

		for (String s : id.split("-"))
		{
			if (s.startsWith(sym + "_")) return s;
		}
		System.out.println("id = " + id);
		System.out.println("sym = " + sym);
		throw new RuntimeException("Should not reach here!");
	}

	private static void readValidationResults() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Aslan-platelet/";
		String valDir = dir + "validation/";
		String blotFile = valDir + "blots.csv";
		String[] header = Files.lines(Paths.get(blotFile)).findFirst().get().split("\t");
		double fdrThr = 0.1;

		Map<String, Set<String>> drugTargetMap = Files.lines(Paths.get(blotFile)).skip(20).limit(8)
			.map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> new HashSet<>(Arrays.asList(t[1].split(", ")))));

		Map<String, Map<String, Boolean>> dirMap = new HashMap<>();
		Files.lines(Paths.get(blotFile)).skip(31).map(l -> l.split("\t")).forEach(t ->
		{
			Map<String, Boolean> map = new HashMap<>();
			double ref = Double.valueOf(t[3]);
			for (int i = 4; i < t.length; i++)
			{
				map.put(header[i], Double.valueOf(t[i]) < ref);
			}
			dirMap.put(t[0], map);
		});


		PhosphoGraph graph = new PhosphoGraph("Validation", SignedType.PHOSPHORYLATES.getTag());

		Files.lines(Paths.get(blotFile)).skip(1).limit(11).map(l -> l.split("\t")).forEach(t ->
		{
			String gene = t[0].substring(0, t[0].indexOf(" "));
			Set<String> sites = new HashSet<>(Arrays.asList(t[0].substring(t[0].indexOf(" ") + 1).split(", ")));

			for (int i = 4; i < t.length; i++)
			{
				if (t[i].startsWith("<")) t[i] = t[i].substring(1);
				double val = Double.valueOf(t[i]);
				if (val <= fdrThr)
				{
					if (!dirMap.get(t[0]).get(header[i]))
					{
						System.out.println("trg = " + t[0] + ", drug = " + header[i]);
					}
					for (String eff : drugTargetMap.get(header[i]))
					{
						graph.putRelation(eff, gene, "", CollectionUtil.merge(sites, ";"));
					}
				}
			}
		});

		graph.write(valDir + "validations.sif");
	}

	private static void compareResultsToValidations() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Aslan-platelet/";
		String valDir = dir + "validation/";
		String resDir = dir + "merged/";

		PhosphoGraph valGraph = new PhosphoGraph("Validation", SignedType.PHOSPHORYLATES.getTag());
		valGraph.load(valDir + "validations.sif", Collections.singleton(valGraph.getEdgeType()));
		PhosphoGraph resGraph = new PhosphoGraph("Results", SignedType.PHOSPHORYLATES.getTag());
		resGraph.load(resDir + "sum-relax1aa.sif", Collections.singleton(resGraph.getEdgeType()));

		Set<String> valSources = valGraph.getOneSideSymbols(true);
		Set<String> valTargets = valGraph.getOneSideSymbols(false);

		for (String source : valSources)
		{
			for (String target : resGraph.getDownstream(source))
			{
				if (valTargets.contains(target))
				{
					System.out.print(source + " ----> " + target + "\t");
					if (valGraph.getDownstream(source).contains(target))
					{
						System.out.println("+");
					}
					else
					{
						System.out.println();
					}
				}
			}
		}
	}
}
