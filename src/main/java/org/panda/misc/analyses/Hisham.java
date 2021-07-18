package org.panda.misc.analyses;

import org.panda.misc.causalpath.CausalPathSubnetwork;
import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.resource.UniProtSequence;
import org.panda.utility.*;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.StouffersCombinedProbability;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

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
//		convertHumanProteomics();
//		checkIntersection();
//		convertMouseProteomics();
//		convertRNAseq();
//		getHumanMouseIntersection();
		generateSubgraphs();

//		addModeratedValues();

//		comparePhosphoPeptidesHumanAndMouse();
	}

	private static void convertMouseProteomics() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-re-redo/";
		String protFile = dir + "PR491_Fusion_Mouse_Proteome.txt";
		String phosphoFile1 = dir + "PR443_QX_EVP_Mouse_PhosphoProteome.txt";
		String phosphoFile2 = dir + "PR607_QX_CT_Mouse_PhosphoProteome.txt";
		String phosphoFile3 = dir + "PR607_QX_EVP_Mouse_PhosphoProteome.txt";
		String header = "ID\tSymbols\tSites\tEffect\tNAIVE-1\tNAIVE-2\tNAIVE-3\tPRIMED-1\tPRIMED-2\tPRIMED-3\tREPRIMED-1\tREPRIMED-2\tREPRIMED-3";

		List<String> outLines = new ArrayList<>();

		readMouseProtFile(protFile, outLines);
		writeConvertedData(dir, outLines, header, "mouse-proteome.tsv");

		outLines.clear();
//		readMousePhosphoFile(phosphoFile1, outLines);
		readMousePhosphoFileUnspecific(phosphoFile1, outLines);
		writeConvertedData(dir, outLines, header, "mouse-phosphoproteome-pr443-evp.tsv");

		outLines.clear();
//		readMousePhosphoFile(phosphoFile2, outLines);
		readMousePhosphoFileUnspecific(phosphoFile2, outLines);
		writeConvertedData(dir, outLines, header, "mouse-phosphoproteome-pr607-ct.tsv");

		outLines.clear();
//		readMousePhosphoFile(phosphoFile3, outLines);
		readMousePhosphoFileUnspecific(phosphoFile3, outLines);
		writeConvertedData(dir, outLines, header, "mouse-phosphoproteome-pr607-evp.tsv");
	}

	private static String chooseOne(Set<String> hSyms, String mSym)
	{
		if (hSyms.size() == 1) return hSyms.iterator().next();

		String cap = mSym.toUpperCase();
		if (hSyms.contains(cap)) return cap;

		return hSyms.iterator().next();
	}

	public static void convertHumanProteomics() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/";
		String protFile = dir + "PR581_Velos_Human_Proteome.txt";
		String phosphoFile1 = dir + "PR596_QX_CT_Human_PhosphoProteome.txt";
		String phosphoFile2 = dir + "PR596_QX_EVP_Human_PhosphoProteome.txt";

		List<String> outLines = new ArrayList<>();
		String header = "ID\tSymbols\tSites\tEffect\tNAIVE-1\tNAIVE-2\tNAIVE-3\tPRIMED-1\tPRIMED-2\tPRIMED-3\tREPRIMED-1\tREPRIMED-2\tREPRIMED-3";

		readHumanProtFile(protFile, outLines);
		writeConvertedData(dir, outLines, header, "human-proteome.tsv");

		outLines.clear();
		readHumanPhosphoFile(phosphoFile1, outLines);
		writeConvertedData(dir, outLines, header, "human-phosphoproteome-ct.tsv");

		outLines.clear();
		readHumanPhosphoFile(phosphoFile2, outLines);
		writeConvertedData(dir, outLines, header, "human-phosphoproteome-evp.tsv");
	}

	private static void writeConvertedData(String dir, List<String> outLines, String header, String s) throws IOException
	{
		keepMostDeviated(outLines);
		Set<String> ids = new HashSet<>();
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + s));
		writer.write(header);
		for (String line : outLines)
		{
			String id = line.substring(0, line.indexOf("\t"));
			if (!ids.contains(id))
			{
				FileUtil.lnwrite(line, writer);
				ids.add(id);
			}
			else
			{
				System.out.println("Duplicate ID skipped: " + id);
			}
		}
		writer.close();
	}

	private static void keepMostDeviated(List<String> outLines)
	{
		Map<String, Double> devMap = new HashMap<>();
		Map<String, String> lineMap = new HashMap<>();

		for (String line : outLines)
		{
			String[] t = line.split("\t");
			String id = t[0];
			double[] d = new double[t.length - 4];

			for (int i = 0; i < d.length; i++)
			{
				d[i] = Math.pow(2, Double.valueOf(t[i + 4]));
			}

			double sd = Summary.stdev(d);

			if (devMap.getOrDefault(id, 0D) < sd)
			{
				devMap.put(id, sd);
				lineMap.put(id, line);
			}
		}

		outLines.retainAll(lineMap.values());
	}

	private static void readHumanPhosphoFile(String phosphoFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(phosphoFile)).findFirst().get().split("\t");
		int[] expInds = new int[9];
		expInds[0] = ArrayUtil.indexOf(header, "NAIVE 01", "NAIVE-01");
		expInds[1] = ArrayUtil.indexOf(header, "NAIVE 02", "NAIVE-02", "NAIVE -02");
		expInds[2] = ArrayUtil.indexOf(header, "NAIVE 03", "NAIVE-03");
		expInds[3] = ArrayUtil.indexOf(header, "PRIMED 01", "PRIMED-01");
		expInds[4] = ArrayUtil.indexOf(header, "PRIMED 02", "PRIMED-02");
		expInds[5] = ArrayUtil.indexOf(header, "PRIMED 03", "PRIMED-03");
		expInds[6] = ArrayUtil.indexOf(header, "ENDODERM 01", "ENDODERM-01");
		expInds[7] = ArrayUtil.indexOf(header, "ENDODERM 02", "ENDODERM-02");
		expInds[8] = ArrayUtil.indexOf(header, "ENDODERM 03", "ENDODERM-03");

		Files.lines(Paths.get(phosphoFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
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

			for (int i : expInds)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			outLines.add(sb.toString());
		});
	}

	private static void readMousePhosphoFile(String phosphoFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(phosphoFile)).findFirst().get().split("\t");
		int[] expInds = new int[9];

		expInds[0] = ArrayUtil.indexOf(header, "zi-rep01", "long term01", "long term 01");
		expInds[1] = ArrayUtil.indexOf(header, "zi-rep02", "long term02", "long term 02");
		expInds[2] = ArrayUtil.indexOf(header, "zi-rep03", "long term03", "long term 03");
		expInds[3] = ArrayUtil.indexOf(header, "serum-rep01", "serum01");
		expInds[4] = ArrayUtil.indexOf(header, "serum-rep02", "serum02");
		expInds[5] = ArrayUtil.indexOf(header, "serum-rep03", "serum03");
		expInds[6] = ArrayUtil.indexOf(header, "EPILC-01", "EPILC 01");
		expInds[7] = ArrayUtil.indexOf(header, "EPILC-02", "EPILC 02");
		expInds[8] = ArrayUtil.indexOf(header, "EPILC-03", "EPILC 03");

		for (int i = 0; i < 6; i++)
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
				sb.append("\t").append(i > 0 ? Math.log(Double.valueOf(t[i])) / Math.log(2) : "NaN");
			}

			outLines.add(sb.toString());
		});
	}

	private static void readMousePhosphoFileUnspecific(String phosphoFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(phosphoFile)).findFirst().get().split("\t");
		int[] expInds = new int[9];

		expInds[0] = ArrayUtil.indexOf(header, "zi-rep01", "long term01", "long term 01");
		expInds[1] = ArrayUtil.indexOf(header, "zi-rep02", "long term02", "long term 02");
		expInds[2] = ArrayUtil.indexOf(header, "zi-rep03", "long term03", "long term 03");
		expInds[3] = ArrayUtil.indexOf(header, "serum-rep01", "serum01");
		expInds[4] = ArrayUtil.indexOf(header, "serum-rep02", "serum02");
		expInds[5] = ArrayUtil.indexOf(header, "serum-rep03", "serum03");
		expInds[6] = ArrayUtil.indexOf(header, "EPILC-01", "EPILC 01");
		expInds[7] = ArrayUtil.indexOf(header, "EPILC-02", "EPILC 02");
		expInds[8] = ArrayUtil.indexOf(header, "EPILC-03", "EPILC 03");

		for (int i = 0; i < 6; i++)
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

			List<String> negSites = sites.stream().map(s -> (s.substring(0, 1) + "-" + s.substring(1))).collect(Collectors.toList());

			Map<String, List<String>> sMap = new HashMap<>();
			for (String hUP : hUPs)
			{
				String symbol = HGNC.get().getSymbol(hUP);
				if (symbol != null)
				{
					sMap.put(symbol, hMap.getOrDefault(hUP, negSites));
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
				sb.append("\t").append(i > 0 ? Math.log(Double.valueOf(t[i])) / Math.log(2) : "NaN");
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

	private static void readHumanProtFile(String protFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(protFile)).findFirst().get().split("\t");
		int[] expInds = new int[9];
		expInds[0] = ArrayUtil.indexOf(header, "NAIVE 01");
		expInds[1] = ArrayUtil.indexOf(header, "NAIVE 02");
		expInds[2] = ArrayUtil.indexOf(header, "NAIVE 03");
		expInds[3] = ArrayUtil.indexOf(header, "PRIMED 01");
		expInds[4] = ArrayUtil.indexOf(header, "PRIMED 02");
		expInds[5] = ArrayUtil.indexOf(header, "PRIMED 03");
		expInds[6] = ArrayUtil.indexOf(header, "ENDODERM 01");
		expInds[7] = ArrayUtil.indexOf(header, "ENDODERM 02");
		expInds[8] = ArrayUtil.indexOf(header, "ENDODERM 03");

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

			for (int i : expInds)
			{
				sb.append("\t").append(Math.log(Double.valueOf(t[i])) / Math.log(2));
			}

			outLines.add(sb.toString());
		});
	}

	private static void readMouseProtFile(String protFile, List<String> outLines) throws IOException
	{
		String[] header = Files.lines(Paths.get(protFile)).findFirst().get().split("\t");
		int[] expInds = new int[9];

		expInds[0] = ArrayUtil.indexOf(header, "zi rep01", "long term 01", "long term-01");
		expInds[1] = ArrayUtil.indexOf(header, "zi rep02", "long term 02", "long term-02");
		expInds[2] = ArrayUtil.indexOf(header, "zi rep03", "long term 03", "long term-03");
		expInds[3] = ArrayUtil.indexOf(header, "serum rep01", "serum 01", "serum-01");
		expInds[4] = ArrayUtil.indexOf(header, "serum rep02", "serum 02", "serum-02");
		expInds[5] = ArrayUtil.indexOf(header, "serum rep03", "serum 03", "serum-03");
		expInds[6] = ArrayUtil.indexOf(header, "EPILC-01");
		expInds[7] = ArrayUtil.indexOf(header, "EPILC-02");
		expInds[8] = ArrayUtil.indexOf(header, "EPILC-03");

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
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/";

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
		String base = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-re-redo/";

//		CausalPathSubnetwork.writeGOIRecursiveCompBased(base + "human", base + "GOI");

		CausalPathSubnetwork.generateNeighborhoodSubgraphsForSignificantsRecursively(base, 0.1);
	}

	private static void checkIntersection() throws IOException
	{
		// human
//		String dir = "/home/ozgun/Analyses/Hisham/Proteome-redo/";
//		Set<String> set1 = Files.lines(Paths.get(dir + "data1.txt")).skip(1).map(l -> l.split("\t")[0]).filter(id -> !id.contains("_")).collect(Collectors.toSet());
//		Set<String> set2 = Files.lines(Paths.get(dir + "data2.txt")).skip(1).map(l -> l.split("\t")[0]).filter(id -> !id.contains("_")).collect(Collectors.toSet());
//
//		CollectionUtil.printNameMapping("Data1", "Data2");
//		CollectionUtil.printVennCounts(set1, set2);

		// mouse
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/";
		Set<String> set1 = Files.lines(Paths.get(dir + "mouse-data1.txt")).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).map(t -> t[0]).collect(Collectors.toSet());
		Set<String> set2 = Files.lines(Paths.get(dir + "mouse-data2.txt")).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).map(t -> t[0]).collect(Collectors.toSet());
		Set<String> set3 = Files.lines(Paths.get(dir + "mouse-data3.txt")).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).map(t -> t[0]).collect(Collectors.toSet());

		CollectionUtil.printNameMapping("Data1", "Data2", "Data3");
		CollectionUtil.printVennCounts(set1, set2, set3);
	}

	private static void comparePhosphoPeptidesHumanAndMouse() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/";

		String hComp = "NAIVE-vs-REPRIMED";
		String mComp = "NAIVE-vs-REPRIMED";

		Map<String, List<Tuple>> human = new HashMap<>();
		Map<String, List<Tuple>> mouse = new HashMap<>();

		load(human, dir + "human-phosphoproteome-evp-moderated.tsv", hComp);
		load(human, dir + "human-phosphoproteome-ct-moderated.tsv", hComp);
		if (!mComp.contains("REPRIMED")) load(mouse, dir + "mouse-phosphoproteome-pr443-evp-moderated.tsv", mComp);
		load(mouse, dir + "mouse-phosphoproteome-pr607-ct-moderated.tsv", mComp);
		load(mouse, dir + "mouse-phosphoproteome-pr607-evp-moderated.tsv", mComp);

		Map<String, Tuple> pvalsHuman = human.keySet().stream().collect(Collectors.toMap(id -> id, id -> combine(human, id)));
		Map<String, Tuple> pvalsMouse = mouse.keySet().stream().collect(Collectors.toMap(id -> id, id -> combine(mouse, id)));

		Map<String, Double> pvals = new HashMap<>();
		pvalsHuman.keySet().forEach(id -> pvals.put(id + "-h", pvalsHuman.get(id).p));
		pvalsMouse.keySet().forEach(id -> pvals.put(id + "-m", pvalsMouse.get(id).p));

		List<String> select = FDR.select(pvals, null, 0.05);

		Set<String> sigHuman = select.stream().filter(id -> id.endsWith("-h")).map(id -> id.substring(0, id.length() - 2)).collect(Collectors.toSet());
		Set<String> sigMouse = select.stream().filter(id -> id.endsWith("-m")).map(id -> id.substring(0, id.length() - 2)).collect(Collectors.toSet());

		Set<String> humanOnly = pvalsHuman.keySet().stream().filter(id -> !pvalsMouse.containsKey(id)).collect(Collectors.toSet());
		Set<String> mouseOnly = pvalsMouse.keySet().stream().filter(id -> !pvalsHuman.containsKey(id)).collect(Collectors.toSet());
		Set<String> common = pvalsMouse.keySet().stream().filter(pvalsHuman::containsKey).collect(Collectors.toSet());

		long cnt_HumOnlyUp = sigHuman.stream().filter(humanOnly::contains).filter(id -> pvalsHuman.get(id).v > 0).count();
		long cnt_HumOnlyDw = sigHuman.stream().filter(humanOnly::contains).filter(id -> pvalsHuman.get(id).v < 0).count();

		long cnt_MouOnlyUp = sigMouse.stream().filter(mouseOnly::contains).filter(id -> pvalsMouse.get(id).v > 0).count();
		long cnt_MouOnlyDw = sigMouse.stream().filter(mouseOnly::contains).filter(id -> pvalsMouse.get(id).v < 0).count();

		long cnt_HumOnlyInsignif = humanOnly.stream().filter(k -> !sigHuman.contains(k)).count();
		long cnt_MouOnlyInsignif = mouseOnly.stream().filter(k -> !sigMouse.contains(k)).count();

		System.out.println("Human mouse compare \"changed\"");
		CollectionUtil.printVennCounts(pvalsHuman.keySet(), pvalsMouse.keySet());
		System.out.println();
		System.out.println("cnt_HumOnlyUp = " + cnt_HumOnlyUp);
		System.out.println("cnt_HumOnlyInsignif = " + cnt_HumOnlyInsignif);
		System.out.println("cnt_HumOnlyDw = " + cnt_HumOnlyDw);
		System.out.println("cnt_MouOnlyUp = " + cnt_MouOnlyUp);
		System.out.println("cnt_MouOnlyInsignif = " + cnt_MouOnlyInsignif);
		System.out.println("cnt_MouOnlyDw = " + cnt_MouOnlyDw);
		System.out.println();

		Set<String> comHumUp = common.stream().filter(sigHuman::contains).filter(id -> pvalsHuman.get(id).v > 0).collect(Collectors.toSet());
		Set<String> comHumDw = common.stream().filter(sigHuman::contains).filter(id -> pvalsHuman.get(id).v < 0).collect(Collectors.toSet());
		Set<String> comMouUp = common.stream().filter(sigMouse::contains).filter(id -> pvalsMouse.get(id).v > 0).collect(Collectors.toSet());
		Set<String> comMouDw = common.stream().filter(sigMouse::contains).filter(id -> pvalsMouse.get(id).v < 0).collect(Collectors.toSet());

		System.out.println("Include up/down distinction for intersecting");
		CollectionUtil.printNameMapping("Human Up", "Human Down", "Mouse Up", "Mouse Down");
		CollectionUtil.printVennCounts(comHumUp, comHumDw, comMouUp, comMouDw);

		Set<String> commonChanged = CollectionUtil.getUnion(comHumUp, comHumDw, comMouUp, comMouDw);

		System.out.println("\ncommonChanged.size() = " + commonChanged.size());
		System.out.println("common unchanged = " + (common.size() - commonChanged.size()));

	}

	private static Tuple combine(Map<String, List<Tuple>> map, String id)
	{
		List<Tuple> list = map.get(id);

		if (list.size() == 1) return list.get(0);

		double[] p = new double[list.size()];
		int[] s = new int[list.size()];

		for (int i = 0; i < p.length; i++)
		{
			Tuple t = list.get(i);
			p[i] = t.p;
			s[i] = t.v >= 0 ? 1 : -1;
		}

		return StouffersCombinedProbability.combineP2Tailed(p, s);
	}

	private static void load(Map<String, List<Tuple>> vals, String filename, String signedPColName) throws IOException
	{
		String[] header = Files.lines(Paths.get(filename)).findFirst().get().split("\t");
		int ind = ArrayUtil.indexOf(header, signedPColName);
		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).forEach(t ->
		{
			String id = t[0].replaceAll("_", "-");

			if (!vals.containsKey(id)) vals.put(id, new ArrayList<>());

			double v = t[ind].startsWith("-") ? -1 : 1;
			double p = Math.abs(Double.valueOf(t[ind]));

			vals.get(id).add(new Tuple(v, p));
		});
	}

	private static void plotPhosphoPeptidesHumanAndMouse() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/";

		String hComp = "NAIVE-vs-REPRIMED";
		String mComp = "NAIVE-vs-REPRIMED";

		Map<String, List<Tuple>> human = new HashMap<>();
		Map<String, List<Tuple>> mouse = new HashMap<>();

		load(human, dir + "human-phosphoproteome-evp-moderated.tsv", hComp);
		load(human, dir + "human-phosphoproteome-ct-moderated.tsv", hComp);
		if (!mComp.contains("REPRIMED")) load(mouse, dir + "mouse-phosphoproteome-pr443-evp-moderated.tsv", mComp);
		load(mouse, dir + "mouse-phosphoproteome-pr607-ct-moderated.tsv", mComp);
		load(mouse, dir + "mouse-phosphoproteome-pr607-evp-moderated.tsv", mComp);

		Map<String, Tuple> pvalsHuman = human.keySet().stream().collect(Collectors.toMap(id -> id, id -> combine(human, id)));
		Map<String, Tuple> pvalsMouse = mouse.keySet().stream().collect(Collectors.toMap(id -> id, id -> combine(mouse, id)));
	}

		private static void addModeratedValues() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-re-redo/";

		Map<String, String> mapD1 = readModeratedValues(dir + "D1.csv");
		Map<String, String> mapD2 = readModeratedValues(dir + "D2.csv");

		String file = dir + "mouse-phosphoproteome-pr443-evp";

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file + "-moderated.tsv"));

		Files.lines(Paths.get(file + ".tsv")).forEach(l ->
		{
			String[] t = l.split("\t");


			if (t[0].equals("ID")) FileUtil.write(l + "\tNAIVE-vs-PRIMED\tNAIVE-vs-REPRIMED", writer);
			else FileUtil.lnwrite(l + "\t" + mapD1.get(t[0]) + "\t" + mapD2.get(t[0]).replaceAll("NA", "NaN"), writer);
		});

		writer.close();
	}

	private static Map<String, String> readModeratedValues(String file) throws IOException
	{
		return FileUtil.readMap(file, ",", "ID", "SignedP");
	}

	private static void convertRNAseq() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Analyses/Hisham/Proteome-redo/rnaseq-data/";

		Map<String, Double> mNVP = FileUtil.lines(dir + "MouseNaiveVsPrime.csv").skip(1).map(l -> l.split(",")).filter(t -> !t[2].equals("NA"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf((t[2].startsWith("-") ? "-" : "") + (t[5].equals("0") ? "1E-200" : t[5]))));

		Map<String, Double> mNVE = FileUtil.lines(dir + "MouseNaiveVsEpiLC.csv").skip(1).map(l -> l.split(",")).filter(t -> !t[2].equals("NA"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf((t[2].startsWith("-") ? "-" : "") + (t[5].equals("0") ? "1E-200" : t[5]))));

		Map<String, Double> hNVP = FileUtil.lines(dir + "HumanNaiveVsPrime.csv").skip(1).map(l -> l.split(",")).filter(t -> !t[5].equals("NA"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf((t[2].startsWith("-") ? "-" : "") + (t[5].equals("0") ? "1E-200" : t[5]))));

		BufferedWriter writerM = Files.newBufferedWriter(Paths.get(dir + "mouse-rnaseq.tsv"));
		BufferedWriter writerH = Files.newBufferedWriter(Paths.get(dir + "human-rnaseq.tsv"));

		writerM.write("ID\tNAIVE-vs-PRIMED\tNAIVE-vs-REPRIMED");
		Set<String> mSyms = new HashSet<>(mNVP.keySet());
		mSyms.addAll(mNVE.keySet());
		mSyms.forEach(mSym ->
		{
			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);
			if (hSyms.isEmpty()) return;

			String hSym = chooseOne(hSyms, mSym);

			FileUtil.lnwrite(hSym + "\t" + mNVP.getOrDefault(mSym, Double.NaN) + "\t" + mNVE.getOrDefault(mSym, Double.NaN), writerM);

		});

		writerH.write("ID\tNAIVE-vs-PRIMED");
		hNVP.keySet().forEach(sym -> FileUtil.lnwrite(sym + "\t" + hNVP.get(sym), writerH));

		writerM.close();
		writerH.close();
	}

	private static Set<String> GOI = new HashSet<>(Arrays.asList((
		"DNMT1\n" +
		"DNMT3A\n" +
		"DNMT3B\n" +
		"DNMT3L\n" +
		"TET1\n" +
		"TET2\n" +
		"TET3\n" +
		"E2F1\n" +
		"EED\n" +
		"SUZ12\n" +
		"EZH1\n" +
		"EZH2\n" +
		"JARID2\n" +
		"KDM1A\n" +
		"WDR5\n" +
		"KMT2A\n" +
		"ASH2L\n" +
		"TRIM28\n" +
		"SETDB1\n" +
		"ZNF809\n" +
		"EHMT2\n" +
		"CBX5\n" +
		"SUV39H1\n" +
		"SUV39H2\n" +
		"ATRX\n" +
		"ZNF274\n" +
		"SETD5\n" +
		"POU5F1\n" +
		"SOX2\n" +
		"MYC\n" +
		"KLF2\n" +
		"KLF4\n" +
		"ZFP42\n" +
		"DPPA3\n" +
		"NANOG\n" +
		"ESRRB\n" +
		"OTX2\n" +
		"NODAL\n" +
		"ZIC3\n" +
		"TCF3\n" +
		"ASNS\n" +
		"PDHA1\n" +
		"CBS\n" +
		"PHGDH\n" +
		"MAT1A\n" +
		"MAT2B\n" +
		"MAT2A\n" +
		"AGA\n" +
		"ASPG\n" +
		"ASRGL1\n" +
		"MTR\n" +
		"MTRR\n" +
		"BHMT\n" +
		"BCAT1\n" +
		"ACAT1\n" +
		"ASS1\n" +
		"FASN\n" +
		"ACSS1\n" +
		"ACSS2\n" +
		"LPIN1\n" +
		"ACACA\n" +
		"ACLY\n" +
		"AGPAT1\n" +
		"EHHADH\n" +
		"ACSL4\n" +
		"ACSL5\n" +
		"ACSL6\n" +
		"ACOT1\n" +
		"ACOT9\n" +
		"HMGCS1\n" +
		"ACAT2\n" +
		"CHSY1\n" +
		"GNPDA1\n" +
		"PGM1\n" +
		"PGM2\n" +
		"PMM1\n" +
		"MAN1C1\n" +
		"HK1\n" +
		"ABCC5\n" +
		"ALDOC\n" +
		"RRM1\n" +
		"RRM2\n" +
		"AK1\n" +
		"TK1\n" +
		"ADARB1\n" +
		"NPR1\n" +
		"NME3\n" +
		"NME4\n" +
		"MAPK1\n" +
		"GSK3B\n" +
		"MTOR\n" +
		"AKT1\n" +
		"AKT2\n" +
		"PRKAA1\n" +
		"PRKAA2").split("\n")));
}
