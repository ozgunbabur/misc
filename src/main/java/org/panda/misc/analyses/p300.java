package org.panda.misc.analyses;

import org.panda.resource.MGI;
import org.panda.resource.PCPathway;
import org.panda.resource.PCPathwayHGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.GeneSetEnrichment;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class p300
{
	public static final String DIR = "/home/ozgun/Documents/Papers/Authoring/p300/";

	public static final String MICROARRAY_FILE = DIR + "microarray.csv";
	public static final String MICROARRAY_HUMAN_FILE = DIR + "tf-analysis/microarray-human.txt";
	public static final String MS_FILE = DIR + "ms.csv";
	public static final String SUMMARY_FILE = DIR + "summary-lists.csv";
	public static final String FORMAT_FILE = DIR + "format.format";
	public static final String ENRICHMENT_FILE = DIR + "pathway-enrichment-protein-and-rna.txt";

	public static void main(String[] args) throws IOException
	{
//		Set<String>[] summaries = readChangeSummaries();
//		Map<String, Double> microarray = readMicroarraySignedPValues();
//		Map<String, Double> msFoldChanges = readMSFoldChanges(summaries);

//		generateFormatFile(microarray, msFoldChanges);
//		convertMicroarrayFileToHuman();

//		Set<String>[] hSets = convertToHuman(summaries);
//		doPathwayEnrichment(Arrays.stream(hSets).flatMap(Collection::stream).collect(Collectors.toSet()));

		plotEnrichment(ENRICHMENT_FILE);

//		checkSomeOverlaps();
	}

	static Map<String, Double> readMicroarraySignedPValues() throws IOException
	{
		return Files.lines(Paths.get(MICROARRAY_FILE)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[2], t -> Math.signum(Double.valueOf(t[16])) * Double.valueOf(t[17])));
	}

	static Set<String>[] readChangeSummaries()
	{
		return new Set[]{
			FileUtil.getTermsInTabDelimitedColumn(SUMMARY_FILE, 1, 0),
			FileUtil.getTermsInTabDelimitedColumn(SUMMARY_FILE, 3, 0),
			FileUtil.getTermsInTabDelimitedColumn(SUMMARY_FILE, 5, 0),
			FileUtil.getTermsInTabDelimitedColumn(SUMMARY_FILE, 7, 0)};
	}

	static Map<String, Double> readMSFoldChanges(Set<String>[] changedSets) throws IOException
	{
		return Files.lines(Paths.get(MS_FILE)).skip(5).map(l -> l.split("\t"))
			.filter(t -> changedSets[0].contains(t[11]) || changedSets[1].contains(t[11]))
			.collect(Collectors.toMap(t -> t[11], t -> (changedSets[0].contains(t[11]) ? 1 : -1) *
				Math.max(Math.max(Double.valueOf(t[4]), Double.valueOf(t[5])), Double.valueOf(t[6]))));
	}

	static void generateFormatFile(Map<String, Double> microarray, Map<String, Double> msFoldChanges) throws IOException
	{
		Map<String, Double> hMicro = convertToHuman(microarray);
		Map<String, Double> hProt = convertToHuman(msFoldChanges);

		ValToColor mCol = new ValToColor(new double[]{-4, 0, 4}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});
		ValToColor pCol = new ValToColor(new double[]{-2, 0, 2}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(FORMAT_FILE));

		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t0 0 0");

		hMicro.keySet().forEach(sym -> FileUtil.lnwrite("node\t" + sym + "\tcolor\t" + mCol.getColorInString(hMicro.get(sym)), writer));
		hProt.keySet().forEach(sym -> FileUtil.lnwrite("node\t" + sym + "\tborderwidth\t2\nnode\t" + sym + "\tbordercolor\t" + pCol.getColorInString(hProt.get(sym)), writer));

		writer.close();
	}

	static double absMax(double v1, double v2)
	{
		return Math.abs(v1) > Math.abs(v2) ? v1 : v2;
	}

	static Map<String, Double> convertToHuman(Map<String, Double> mMap)
	{
		Map<String, Double> hMap = new HashMap<>();

		for (String mSym : mMap.keySet())
		{
			for (String hSym : MGI.get().getCorrespondingHumanSymbols(mSym))
			{
				if (hMap.containsKey(hSym)) hMap.put(hSym, absMax(hMap.get(hSym), mMap.get(mSym)));
				else hMap.put(hSym, mMap.get(mSym));
			}
		}
		return hMap;
	}

	static Set<String>[] convertToHuman(Set<String>[] summary)
	{
		Set<String>[] sets = new Set[summary.length];

		for (int i = 0; i < summary.length; i++)
		{
			sets[i] = summary[i].stream().map(MGI.get()::getCorrespondingHumanSymbols).flatMap(Collection::stream)
				.collect(Collectors.toSet());
		}

		CollectionUtil.removeIntersection(sets[0], sets[1]);
		CollectionUtil.removeIntersection(sets[2], sets[3]);

		return sets;
	}

	static void convertMicroarrayFileToHuman() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(MICROARRAY_HUMAN_FILE));
		writer.write("ID\tWT_1\tWT_2\tWT_3\tWT_4\tKO_1\tKO_2\tKO_3\tKO_4");
		Files.lines(Paths.get(MICROARRAY_FILE)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String mSym = t[2];
			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);
			if (hSyms.size() == 1)
			{
				String hSym = hSyms.iterator().next();

				FileUtil.lnwrite(hSym + "\t" + ArrayUtil.getString("\t", log(t[5]), log(t[6]), log(t[7]),
					log(t[8]), log(t[9]), log(t[10]), log(t[11]), log(t[12])), writer);
			}
		});

		writer.close();
	}

	static void doPathwayEnrichment(Set<String> genes)
	{
		PCPathwayHGNC pcp = new PCPathwayHGNC();
		Map<String, Double> pvals = pcp.calculateEnrichment(genes, 5, 200);
		Map<String, Double> qvals = FDR.getQVals(pvals, null);
		pvals.keySet().stream().sorted(Comparator.comparing(pvals::get)).filter(k -> pvals.get(k) < 0.05).forEach(k ->
			System.out.println(k + "\t" + pcp.getPathwayName(k) + "\t" + pvals.get(k) + "\t" + qvals.get(k) + "\t" +
				pcp.getOverlappingGenes(k, genes)));
	}


	static void plotEnrichment(String file) throws IOException
	{
		Set<String>[] summaries = readChangeSummaries();
		Set<String>[] hSets = convertToHuman(summaries);
		Set<String> upSet = new HashSet<>(hSets[0]);
		upSet.addAll(hSets[2]);
		Set<String> dwSet = new HashSet<>(hSets[1]);
		dwSet.addAll(hSets[3]);

		String fileWOExt = file.substring(0, file.lastIndexOf("."));
		Map<String, Double> pvals = new HashMap<>();
		Map<String, Set<String>> contributers = new HashMap<>();

		Files.lines(Paths.get(file)).map(l -> l.split("\t"))/*.filter(t -> Double.valueOf(t[2]) < 1)*/.forEach(t ->
		{
			pvals.put(t[1], Double.valueOf(t[2]));
			contributers.put(t[1], new HashSet<>(Arrays.asList(t[4].substring(1, t[4].length() - 1).split(", "))));
		});

		GeneSetEnrichment.graphGeneSetOverlaps(pvals, contributers, fileWOExt, upSet, dwSet);
	}

	static void checkSomeOverlaps()
	{
		Set<String>[] summaries = readChangeSummaries();
		Set<String>[] hSets = convertToHuman(summaries);

//		System.out.println(PCPathway.get().getName("http://identifiers.org/reactome/R-HSA-72706"));

		Set<String> query = PCPathway.get().getGenes("http://identifiers.org/reactome/R-HSA-156902");
//		Set<String> query = new HashSet<>(Arrays.asList("EIF4A2, RPL4, EIF4A1, RPL5, RPL30, RPL3, RPL32, RPL31, RPL34, APEH, RPL8, RPL10A, RPL9, RPL6, RPL7, RPS15, EEF1B2, RPS14, RPS17, RPS16, RPL18A, RPS19, RPS18, RPL36, RPL35, RPL38, RPL37, RPS11, RPL39, RPS10, RPS13, RPS12, RPS9, RPS7, RPL21, EIF1AX, RPS8, RPL23, RPS5, RPS6, RPL22, RPSA, RPL9P7, RPL9P8, RPL9P9, EEF1A1, RPL24, RPL27, RPL26, RPL29, UBA52, RPL28, CHEBI:33715, RPS4Y1, GSPT2, EIF2B5, EIF2B4, RPL41, EIF2B3, EIF2B2, EIF2S2, EEF2, EIF2S1, RPL3L, RPS26, RPS25, RPS28, RPS27, EIF5, EIF2S3, RPS29, RPL27A, CHEBI:33695, RPS20, RPS21, RPS24, EIF4G1, RPS23, RPLP1, RPLP0, CHEBI:45441, RPS4X, RPL7A, RPLP2, EIF2B1, EIF5B, RPL13A, RPS3A, CHEBI:13677, EEF1G, CHEBI:83690, EEF1D, RPL37A, ETF1, RPL10, RPL12, RPL11, RPL36A, CHEBI:14911, RPS15A, CHEBI:16541, EIF4H, RPS3, RPL14, RPL13, RPL15, RPS2, RPL18, RPS27A, EIF4E, RPL17, EIF4B, RPL19, RPL35A, RPL23A, EIF3M, EIF3K, EIF3L, EIF3I, EIF3J, EIF3G, CHEBI:36080, EIF3H, EIF3E, FAU, EIF3F, EIF3C, EIF3D, RPL26L1, EIF3A, EIF3B".split(", ")));
		System.out.println("query.size() = " + query.size());

		for (Set<String> hSet : hSets)
		{
			List<String> intersection = CollectionUtil.getIntersection(hSet, query).stream().sorted().collect(Collectors.toList());
			System.out.println("intersection = " + intersection);
		}
	}

	static String log(String s)
	{
		return String.valueOf(Math.log(Double.valueOf(s)) / Math.log(2));
	}
}
