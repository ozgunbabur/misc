package org.panda.misc.analyses;

import org.panda.causalpath.run.CausalPath;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Histogram;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PNNLProteomics
{
	public static final String DIR = "/home/babur/Documents/TCGA/OV/";
	public static final String OUT_DIR = "/home/babur/Documents/RPPA/TCGA/PNNL/";
	public static final String PROT_FILE = DIR + "PNNL-proteome.txt";
	public static final String PHOSPHO_FILE = DIR + "PNNL-phosphoproteome.txt";
	public static final String OUTPUT_FILE = OUT_DIR + "PNLL-causality-formatted.txt";
	public static final String SUBTYPE_DIR = OUT_DIR + "/subtypes";
	public static final String CLINICAL_FILE = OUT_DIR + "/clinical-data.txt";
	public static final String RNA_SUBTYPES_FILE = SUBTYPE_DIR + "/rna-subtype-map.txt";
	public static final String PROTEOMICS_SUBTYPES_FILE = SUBTYPE_DIR + "/proteomic-subtype-map.txt";

	private final static String COMPARISON_PARAMETERS_PREFIX = "proteomics-platform-file = ../../PNLL-causality-formatted.txt\n" +
		"proteomics-values-file = ../../PNLL-causality-formatted.txt\n" +
		"id-column = ID\n" +
		"symbols-column = Symbols\n" +
		"sites-column = Sites\n" +
		"effect-column = Effect\n" +
		"do-log-transform = false\n" +
		"\n" +
		"threshold-for-data-significance = 0.05\n" +
		"value-transformation = significant-change-of-mean\n" +
		"\n" +
		"relation-filter-type = phospho-primary-expression-secondary\n" +
		"\n" +
		"calculate-network-significance = true\n" +
		"permutations-for-significance = 10000\n" +
		"\n" +
		"tcga-directory = ../..\n" +
		"mutation-effect-file = ../../mutation-effects.txt";


	public static void main(String[] args) throws IOException
	{
//		prepareDataFile();
		prepareSubtypeFolders(RNA_SUBTYPES_FILE, "Subtype-");
//		prepareSubtypeFolders(PROTEOMICS_SUBTYPES_FILE, "Proteomic-Subtype-");
//		preparePlainumStatus();
//		prepareOldVersusYoung();
//		prepareShortVsLongSurvival();
	}



	public static void prepareShortVsLongSurvival() throws IOException
	{
		Map<String, Integer> map = Files.lines(Paths.get(CLINICAL_FILE)).skip(1).map(l -> l.split("\t"))
			.filter(t -> !t[10].startsWith("N")).collect(Collectors.toMap(t -> t[0], t -> Integer.valueOf(t[10])));

		Histogram h = new Histogram(100);
		h.setBorderAtZero(true);
		map.values().stream().forEach(v -> h.count(v));
		h.print();

		String subDir = SUBTYPE_DIR + File.separator + "ShortVsLongSurvival";
		Files.createDirectories(Paths.get(subDir));
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(subDir + File.separator + "parameters.txt"));

		FileUtil.write(COMPARISON_PARAMETERS_PREFIX, writer);

		map.keySet().stream().forEach(s ->
			FileUtil.lnwrite((map.get(s) <= 1000 ? "test-value-column = " : "control-value-column = ") +
				s, writer));

		writer.close();
		CausalPath.main(new String[]{subDir});
	}

	public static void prepareOldVersusYoung() throws IOException
	{
		Map<String, Integer> map = Files.lines(Paths.get(CLINICAL_FILE)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Integer.valueOf(t[3])));

		Histogram h = new Histogram(10);
		h.setBorderAtZero(true);
		map.values().stream().forEach(v -> h.count(v));
		h.print();

		String subDir = SUBTYPE_DIR + File.separator + "OldVsYoung";
		Files.createDirectories(Paths.get(subDir));
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(subDir + File.separator + "parameters.txt"));

		FileUtil.write(COMPARISON_PARAMETERS_PREFIX, writer);

		map.keySet().stream().forEach(s ->
			FileUtil.lnwrite((map.get(s) >= 50 ? "test-value-column = " : "control-value-column = ") +
				s, writer));

		writer.close();
		CausalPath.main(new String[]{subDir});
	}

	public static void preparePlainumStatus() throws IOException
	{
		Map<String, String> map = Files.lines(Paths.get(CLINICAL_FILE)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> t[11]));

		String subDir = SUBTYPE_DIR + File.separator + "PlatinumResistantVsSensitive";
		Files.createDirectories(Paths.get(subDir));
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(subDir + File.separator + "parameters.txt"));

		FileUtil.write(COMPARISON_PARAMETERS_PREFIX, writer);

		map.keySet().stream().filter(s -> !map.get(s).equals("Not available")).forEach(s ->
			FileUtil.lnwrite((map.get(s).equals("Resistant") ? "test-value-column = " : "control-value-column = ") +
				s, writer));

		writer.close();
		CausalPath.main(new String[]{subDir});
	}

	public static void prepareSubtypeFolders(String SUBTYPES_FILE, String prefix) throws IOException
	{
		System.out.println("prefix = " + prefix);
		Set<String> samplesInPNNL = Arrays.stream(Files.lines(Paths.get(OUTPUT_FILE)).findFirst().get().split("\t"))
			.skip(4).collect(Collectors.toSet());

		Map<String, String> idToType = Files.lines(Paths.get(SUBTYPES_FILE)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> t[1]));

		Set<String> types = idToType.values().stream().distinct().collect(Collectors.toSet());

		int i = 0;
		for (String type : types)
		{
//			if (i++ % 4 != 3) continue;

			String subDir = SUBTYPE_DIR + File.separator + prefix + type;
			Files.createDirectories(Paths.get(subDir));
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(subDir + File.separator + "parameters.txt"));

			FileUtil.write(COMPARISON_PARAMETERS_PREFIX, writer);

			samplesInPNNL.stream().filter(s -> idToType.keySet().contains(s)).forEach(s ->
				FileUtil.lnwrite((idToType.get(s).equals(type) ? "test-value-column = " : "control-value-column = ") +
				s, writer));

			writer.close();

			System.out.println("\ntype = " + type);
			CausalPath.main(new String[]{subDir});
		}
	}

	public static void prepareDataFile() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(OUTPUT_FILE));
		writer.write("ID\tSymbols\tSites\tEffect");

		// read total protein

		String[] header = Files.lines(Paths.get(PROT_FILE)).findFirst().get().split("\t");
		Map<String, Set<Integer>> idToIndices = getSampleToIndicesMap(header);
		List<String> ids = new ArrayList<>(idToIndices.keySet());
		Collections.sort(ids);
		ids.forEach(id -> FileUtil.write("\t" + id, writer));

		Files.lines(Paths.get(PROT_FILE)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (t[1].isEmpty()) return;

			FileUtil.lnwrite(t[1] + "\t" + t[1] + "\t\t", writer);

			for (String id : ids)
			{
				FileUtil.write("\t" + getAveragedValue(t, idToIndices.get(id)), writer);
			}
		});

		for (String id : ids)
		{
			System.out.println("value-column = " + id);
		}

		// read phosphoprotein

		header = Files.lines(Paths.get(PHOSPHO_FILE)).findFirst().get().split("\t");
		Map<String, Integer> idToIndex = new HashMap<>();
		for (int i = 5; i < header.length; i++)
		{
			idToIndex.put(header[i], i);
		}

		Set<String> missing = new HashSet<>(idToIndex.keySet());
		missing.removeAll(idToIndices.keySet());
		System.out.println("Missing from total proteins = " + missing);

		Map<String, Map<String, Set<String[]>>> data = new HashMap<>();

		Files.lines(Paths.get(PHOSPHO_FILE)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[0];
			if (sym.isEmpty()) return;
			String site = t[3];
			if (site.isEmpty()) return;
			if (!data.containsKey(sym)) data.put(sym, new HashMap<>());
			if (!data.get(sym).containsKey(site)) data.get(sym).put(site, new HashSet<>());

			data.get(sym).get(site).add(t);
		});

		data.keySet().stream().sorted().forEach(sym ->
			data.get(sym).keySet().stream().sorted(
				(s1, s2) -> Integer.valueOf(s1.substring(s1.lastIndexOf("-") + 2, s1.length() - 1).split("s|t|y")[0])
				.compareTo(Integer.valueOf(s2.substring(s2.lastIndexOf("-") + 2, s2.length() - 1).split("s|t|y")[0]))).forEach(site ->
			{
				FileUtil.lnwrite(site + "\t" + sym + "\t", writer);
				String[] s = site.substring(site.lastIndexOf("-") + 1, site.length() - 1).split("s|t|y");
				FileUtil.write(ArrayUtil.getString("|", s) + "\t", writer);
				Set<String[]> rows = data.get(sym).get(site);

				ids.forEach(id ->
				{
					if (idToIndex.containsKey(id))
					{
						int i = idToIndex.get(id);
						Set<String> vals = rows.stream().map(r -> r.length > i ? r[i] : "").collect(Collectors.toSet());
						FileUtil.write("\t" + getMax(vals), writer);
					}
					else FileUtil.write("\tNaN", writer);
				});
			}));

		writer.close();
	}

	private static Map<String, Set<Integer>> getSampleToIndicesMap(String[] header)
	{
		Map<String, Set<Integer>> idToIndex = new HashMap<>();
		for (int i = 2; i < header.length; i++)
		{
			assert header[i].startsWith("PNNL-") || header[i].startsWith("JHU-");

			String id = header[i].substring(header[i].indexOf("-") + 1);

			if (!idToIndex.containsKey(id)) idToIndex.put(id, new HashSet<>());
			idToIndex.get(id).add(i);
		}
		return idToIndex;
	}

	private static double getAveragedValue(String[] row, Set<Integer> indices)
	{
		List<Double> vals = new ArrayList<>();
		for (Integer index : indices)
		{
			if (row.length > index && !row[index].isEmpty()) vals.add(Double.valueOf(row[index]));
		}
		if (vals.isEmpty()) return Double.NaN;
		else return ArrayUtil.mean(vals.toArray(new Double[vals.size()]));
	}

	private static double getMax(Set<String> vals)
	{
		Optional<Double> val = vals.stream().filter(s -> !s.isEmpty()).map(Double::valueOf).max(Double::compare);
		if (val.isPresent()) return val.get();
		else return Double.NaN;
	}
}
