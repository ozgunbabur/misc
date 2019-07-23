package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.resource.UniProtSequence;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.graph.UndirectedGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Hisham
{
	public static void main(String[] args) throws IOException
	{
//		convertProteomics();
//		checkIntersection();
//		convertTranscriptomics();
//		convertMouseProteomics();
//		convertMouseTranscriptomics();
		getHumanMouseIntersection();
//		generateSubgraphs();
	}

	private static void convertMouseProteomics() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";
		String protFile1 = dir + "PR443_QX_Mouse_Proteome.txt";
		String protFile2 = dir + "PR491_Fusion_Mouse_Proteome.txt";
		String protFile3 = dir + "PR491_QX_Mouse_Proteome.txt";
		String phosphoFile1 = dir + "PR443_QX_EVP_Mouse_PhosphoProteome.txt";
		String phosphoFile2 = dir + "PR607_QX_CT_Mouse_PhosphoProteome.txt";
		String phosphoFile3 = dir + "PR607_QX_EVP_Mouse_PhosphoProteome.txt";

		List<String> outLines = new ArrayList<>();

		readMousePhosphoFile(phosphoFile3, outLines);
		readMouseProtFile(protFile3, outLines);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "mouse-data3.txt"));
		writer.write("ID\tSymbols\tSites\tEffect\tNAIVE-01\tNAIVE-02\tNAIVE-03\tPRIMED-01\tPRIMED-02\tPRIMED-03");

		outLines.forEach(l -> FileUtil.lnwrite(l, writer));

		writer.close();
	}

	private static void convertTranscriptomics() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";
		String inFile = dir + "human_transcriptome_data.txt";
		String outFile = dir + "rnaseq.txt";

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		writer.write("Symbol\tNAIVE-01\tNAIVE-02\tNAIVE-03\tPRIMED-01\tPRIMED-02\tPRIMED-03");

		Files.lines(Paths.get(inFile)).skip(1).filter(l -> !l.isEmpty()).map(l -> l.split("\t")).forEach(t ->
			FileUtil.lnwrite(t[0] + "\t" + t[12] + "\t" + t[13] + "\t" + t[14] + "\t" + t[15] + "\t" + t[16] +
				"\t" + t[17], writer));

		writer.close();
	}

	private static void convertMouseTranscriptomics() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";
		String inFile = dir + "E14Ser2i_RNASeq_RPKM.txt";
		String outFile = dir + "mouse-rnaseq.txt";

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		writer.write("Symbol\tNAIVE-01\tNAIVE-02\tNAIVE-03\tPRIMED-01\tPRIMED-02\tPRIMED-03");

		Files.lines(Paths.get(inFile)).skip(1).filter(l -> !l.isEmpty()).map(l -> l.split("\t")).forEach(t ->
		{
			String mSym = t[0];

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);
			if (hSyms.isEmpty()) return;

			String hSym = chooseOne(hSyms, mSym);

			FileUtil.lnwrite(hSym + "\t" + t[12] + "\t" + t[13] + "\t" + t[14] + "\t" + t[18] + "\t" + t[19] +
				"\t" + t[20], writer);
		});

		writer.close();
	}

	private static String chooseOne(Set<String> hSyms, String mSym)
	{
		if (hSyms.size() == 1) return hSyms.iterator().next();

		String cap = mSym.toUpperCase();
		if (hSyms.contains(cap)) return cap;

		return hSyms.iterator().next();
	}

	public static void convertProteomics() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";
		String protFile1 = dir + "PR581_QX_Human_Proteome.txt";
		String protFile2 = dir + "PR581_Velos_Human_Proteome.txt";
		String phosphoFile1 = dir + "PR596_QX_CT_Human_PhosphoProteome.txt";
		String phosphoFile2 = dir + "PR596_QX_EVP_Human_PhosphoProteome.txt";

		List<String> outLines = new ArrayList<>();

		readPhosphoFile(phosphoFile2, outLines);
		readProtFile(protFile2, outLines);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "data2.txt"));
		writer.write("ID\tSymbols\tSites\tEffect");
		String[] t = Files.lines(Paths.get(phosphoFile1)).findFirst().get().split("\t");
		for (int i = t.length - 9; i < t.length; i++)
		{
			writer.write("\t" + t[i]);
		}

		outLines.forEach(l -> FileUtil.lnwrite(l, writer));

		writer.close();
	}

	private static void readPhosphoFile(String phosphoFile1, List<String> outLines) throws IOException
	{
		Files.lines(Paths.get(phosphoFile1)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String seq = t[1].substring(4, t[1].length() - 4);
			seq = seq.replaceAll("\\.", "");

			List<Integer> localPos = getLocalPhosphoLocs(t[2]);

			int startLoc = UniProtSequence.get().getStartLocation(t[5], seq);
			if (startLoc < 1)
			{
				System.out.println("Cannot match peptide: " + t[5] + "\t" + seq);
				return;
			}

			List<String> sites = getGlobalSites(seq, localPos, startLoc);
			if (sites == null) return;

			String symbol = HGNC.get().getSymbol(t[5]);

			if (symbol == null)
			{
				System.out.println("Cannot find symbol of \"" + t[5] + "\".");
				return;
			}

			StringBuilder sb = new StringBuilder();

			String id = symbol;
			for (String site : sites)
			{
				id += "_" + site;
			}
			sb.append(id).append("\t").append(symbol).append("\t").append(CollectionUtil.merge(sites, "|")).append("\t");

			for (int i = 7; i < t.length; i++)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			outLines.add(sb.toString());
		});
	}

	private static void readMousePhosphoFile(String phosphoFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(phosphoFile)).findFirst().get().split("\t");
		int[] expInds = new int[6];

		expInds[0] = ArrayUtil.indexOf(header, "serum-rep01", "serum01");
		expInds[1] = ArrayUtil.indexOf(header, "serum-rep02", "serum02");
		expInds[2] = ArrayUtil.indexOf(header, "serum-rep03", "serum03");
		expInds[3] = ArrayUtil.indexOf(header, "zi-rep01", "long term01", "long term 01");
		expInds[4] = ArrayUtil.indexOf(header, "zi-rep02", "long term02", "long term 02");
		expInds[5] = ArrayUtil.indexOf(header, "zi-rep03", "long term03", "long term 03");

		for (int i = 0; i < expInds.length; i++)
		{
			if (expInds[i] < 0)
			{
				System.err.println("Cannot find a column: " + i);
				return;
			}
		}

		Files.lines(Paths.get(phosphoFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String seq = t[1].substring(t[1].indexOf(".") + 1, t[1].lastIndexOf("."));

			List<Integer> localPos = getLocalPhosphoLocs(t[2]);

			String mID = t[5];

			int startLoc = UniProtSequence.get().getStartLocation(mID, seq);
			if (startLoc < 1)
			{
				System.out.println("Cannot match peptide: " + mID + "\t" + seq);
				return;
			}

			List<String> sites = getGlobalSites(seq, localPos, startLoc);
			if (sites == null) return;

			StringBuilder sb = new StringBuilder();

			String mSym = MGI.get().getSymbol(mID);
			if (mSym == null) return;

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);
			Set<String> hUPs = hSyms.stream().map(HGNC.get()::getUniProt).filter(Objects::nonNull).collect(Collectors.toSet());

			Map<String, List<String>> hMap = SiteMappingMouseToHuman.get().mapToHumanSite(mID, sites.toArray(new String[sites.size()]));

			new HashSet<>(hMap.keySet()).forEach(up ->
			{
				if (!hUPs.contains(up)) hMap.remove(up);
			});

			if (hMap.isEmpty())
			{
				System.out.println("Cannot find corresponding human prot: " + mID);
				return;
			}

			Map<String, List<String>> sMap = new HashMap<>();
			for (String hUP : hMap.keySet())
			{
				String symbol = HGNC.get().getSymbol(hUP);
				if (symbol != null)
				{
					sMap.put(symbol, hMap.get(hUP));
				}
			}

			if (sMap.isEmpty())
			{
				System.out.println("Cannot find symbol of human protein(s): " + hMap.keySet());
				return;
			}

			String syms = CollectionUtil.merge(sMap.keySet().stream().sorted().collect(Collectors.toList()), " ");
			String siteStr = CollectionUtil.merge(sMap.keySet().stream().sorted().map(k -> CollectionUtil.merge(sMap.get(k), "|")).collect(Collectors.toList()), " ");
			String id = CollectionUtil.merge(sMap.keySet().stream().sorted().map(k -> k + "-" + CollectionUtil.merge(sMap.get(k), "-")).collect(Collectors.toList()), "_");

			sb.append(id).append("\t").append(syms).append("\t").append(siteStr).append("\t");

			for (int i : expInds)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			outLines.add(sb.toString());
		});
	}

	private static List<String> getGlobalSites(String seq, List<Integer> localPos, int startLoc)
	{
		List<String> sites = new ArrayList<>();
		for (Integer pos : localPos)
		{
			String aa = seq.substring(pos - 1, pos);

			if (!aa.equals("S") && !aa.equals("T") && !aa.equals("Y"))
			{
				System.out.println("aa = " + aa);
			}

			sites.add(aa + (startLoc + pos - 1));
		}

		if (sites.isEmpty()) return null;
		return sites;
	}

	private static void readProtFile(String protFile, List<String> outLines) throws IOException
	{
		Files.lines(Paths.get(protFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String symbol = HGNC.get().getSymbol(t[2]);

			if (symbol == null)
			{
				if (t[3].contains("GN="))
				{
					symbol = t[3].substring(t[3].indexOf("GN=") + 3, t[3].indexOf(" ", t[3].indexOf("GN=")));
				}

				if (symbol == null)
				{
					System.out.println("Cannot find symbol of \"" + t[2] + "\". File has \"" + t[3] + "\"");
					return;
				}
			}

			StringBuilder sb = new StringBuilder();

			sb.append(symbol).append("\t").append(symbol).append("\t\t");

			// reorganizing columns to match phospho files

			for (int i = 11; i < 14; i++)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			for (int i = 8; i < 11; i++)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			for (int i = 14; i < t.length; i++)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			outLines.add(sb.toString());
		});
	}

	private static void readMouseProtFile(String protFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(protFile)).findFirst().get().split("\t");
		int[] expInds = new int[6];

		expInds[0] = ArrayUtil.indexOf(header, "serum rep01", "serum 01", "serum-01");
		expInds[1] = ArrayUtil.indexOf(header, "serum rep02", "serum 02", "serum-02");
		expInds[2] = ArrayUtil.indexOf(header, "serum rep03", "serum 03", "serum-03");
		expInds[3] = ArrayUtil.indexOf(header, "zi rep01", "long term 01", "long term-01");
		expInds[4] = ArrayUtil.indexOf(header, "zi rep02", "long term 02", "long term-02");
		expInds[5] = ArrayUtil.indexOf(header, "zi rep03", "long term 03", "long term-03");

		for (int i = 0; i < expInds.length; i++)
		{
			if (expInds[i] < 0)
			{
				System.err.println("Cannot find a column: " + i);
				return;
			}
		}

		Files.lines(Paths.get(protFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String mUP = t[2];

			String mSym = MGI.get().getSymbol(mUP);
			if (mSym == null) return;

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);

			if (hSyms.isEmpty())
			{
				System.out.println("Cannot find human symbol for: " + mUP);
				return;
			}

			StringBuilder sb = new StringBuilder();

			sb.append(CollectionUtil.merge(hSyms, "_")).append("\t").append(CollectionUtil.merge(hSyms, " ")).append("\t\t");

			// write experiment readouts

			for (int i : expInds)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			outLines.add(sb.toString());
		});
	}

	private static List<Integer> getLocalPhosphoLocs(String str)
	{
		List<Integer> list = new ArrayList<>();

		int i = str.indexOf("Phospho [");

		do
		{
			int j = str.indexOf("]", i + 9);

			String s = str.substring(i + 9, j);

			for (String t : s.split("; "))
			{
				if (t.contains("("))
				{
					Integer loc = Integer.valueOf(t.substring(1, t.indexOf("(")));
					Double score = Double.valueOf(t.substring(t.indexOf("(") + 1, t.indexOf(")")));

					if (score >= 20) list.add(loc);
					else System.out.println("Low prob phospho = " + str);
				}
			}

			i = str.indexOf("Phospho [", i + 9);
		}
		while(i > 0);

		return list;
	}

	private static void getHumanMouseIntersection() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";

		String humDir, mouDir;

		humDir = dir + "primed-vs-naive-phospho/";
		mouDir = dir + "mouse-primed-vs-naive-phospho/";
		SIFFileUtil.writeIntersection(humDir + "causative.sif", mouDir + "causative.sif", humDir + "causative-mouse-intersect.sif");

		humDir = dir + "primed-vs-naive-expression/";
		mouDir = dir + "mouse-primed-vs-naive-expression/";
		SIFFileUtil.writeIntersection(humDir + "causative.sif", mouDir + "causative.sif", humDir + "causative-mouse-intersect.sif");

		humDir = dir + "primed-vs-naive-expression-rnaseq/";
		mouDir = dir + "mouse-primed-vs-naive-expression-rnaseq/";
		SIFFileUtil.writeIntersection(humDir + "causative.sif", mouDir + "causative.sif", humDir + "causative-mouse-intersect.sif");
	}


	private static void generateSubgraphs() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";

		String humDir, mouDir;

		humDir = dir + "primed-vs-naive-expression-rnaseq/";

		Set<String> seed = new HashSet<>(Arrays.asList("KLF17", "NANOG", "KLF2", "KLF4", "DPPA2", "PRDM14", "GATA6", "DNMT3L"));

		seed.stream().sorted().forEach(s -> System.out.print(", " + s));

		UndirectedGraph graph = (UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.INTERACTS_WITH);
		graph.merge((UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH));


		SIFFileUtil.generateSubsetAndAddPPI(humDir + "causative.sif", seed, humDir + "causative-subset.sif", graph);
	}

	private static void checkIntersection() throws IOException
	{
		// human
//		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";
//		Set<String> set1 = Files.lines(Paths.get(dir + "data1.txt")).skip(1).map(l -> l.split("\t")[0]).filter(id -> !id.contains("_")).collect(Collectors.toSet());
//		Set<String> set2 = Files.lines(Paths.get(dir + "data2.txt")).skip(1).map(l -> l.split("\t")[0]).filter(id -> !id.contains("_")).collect(Collectors.toSet());
//
//		CollectionUtil.printNameMapping("Data1", "Data2");
//		CollectionUtil.printVennCounts(set1, set2);

		// mouse
		String dir = "/home/ozgun/Analyses/Hisham/Proteome/";
		Set<String> set1 = Files.lines(Paths.get(dir + "mouse-data1.txt")).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).map(t -> t[0]).collect(Collectors.toSet());
		Set<String> set2 = Files.lines(Paths.get(dir + "mouse-data2.txt")).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).map(t -> t[0]).collect(Collectors.toSet());
		Set<String> set3 = Files.lines(Paths.get(dir + "mouse-data3.txt")).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).map(t -> t[0]).collect(Collectors.toSet());

		CollectionUtil.printNameMapping("Data1", "Data2", "Data3");
		CollectionUtil.printVennCounts(set1, set2, set3);
	}
}
