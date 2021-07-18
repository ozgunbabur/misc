package org.panda.misc.analyses;

import org.panda.resource.HGNC;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public class Hasmik
{
	public static final String DATA_DIR = "/Users/ozgun/Documents/Data/Hasmik/";
	public static final String BASE = "/Users/ozgun/Documents/Analyses/SigPath-CausalPath/";

	public static void main(String[] args) throws IOException
	{
//		convertPlatform();
//		convertMedulo();
//		convertMeduloNoNorm();
//		convertPDX();
		convertPDXDrugT();
//		convertPerturbation();
//		prepareAnalysisFolders();

//		printMissingSites(BASE + "Medulloblastoma");
//		System.out.println("-----");
//		printMissingSites(BASE + "PDX");
//		System.out.println("-----");
//		printMissingSites(BASE + "Perturbation");
	}

	private static void convertPlatform() throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(BASE + "platform.txt");
		writer.write("ID\tSymbols\tSites\tEffect\tNotes");
		Set<String> idMem = new HashSet<>();
		String rawFile = DATA_DIR + "platform-raw.csv";
		String[] header = FileUtil.readHeader(rawFile, 3);
		int idIndex = ArrayUtil.indexOf(header, "pSite ID for Causal Path analysis");
		int symIndex = ArrayUtil.indexOf(header, "Gene symbol");
		int siteIndex = ArrayUtil.indexOf(header, "phosphosite");
		int noteIndex = ArrayUtil.indexOf(header, "Rational for inclusion  of nominated sites");

		FileUtil.lines(rawFile).skip(4).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[idIndex].trim();
			if (idMem.contains(id)) return;
			idMem.add(id);

			String syms = ArrayUtil.merge(" ", t[symIndex].split("/"));

			for (String sym : syms.split(" "))
			{
				if (!sym.equals(HGNC.get().getSymbol(sym)))
				{
					System.out.println("sym to fix = " + sym);
				}
			}

			String sites = ArrayUtil.merge(" ", t[siteIndex].replaceAll("_", "|").split("/"));
			if (sites.endsWith("y ")) sites = sites.substring(0, sites.length() - 2);

			String notes = t[noteIndex];
			String effect =
				notes.toLowerCase().contains("activation mark") ||
				notes.toLowerCase().contains("mark of activation") ||
					notes.toLowerCase().contains("increased kinase activity")	? "a" :
				notes.toLowerCase().contains("inhibitory site") ||
					notes.toLowerCase().contains("mark of inactivation") ? "i" : "";

			FileUtil.lnwrite(ArrayUtil.merge("\t", id, syms, sites, effect, notes), writer);
		});
		writer.close();
	}

	private static void convertMedulo() throws IOException
	{
		String inFile = "medullo.csv";
		int idIndex = 0;
		String outFile = "data-medullo.txt";
		int[] pInds = new int[]{1, 2, 3};
		int[] fcInds = new int[]{10, 11, 12};

		convert(inFile, outFile, idIndex, pInds, fcInds);
	}

	private static void convertMeduloNoNorm() throws IOException
	{
		String inFile = "medullo-no-norm.csv";
		int idIndex = 0;
		String outFile = "data-medullo-no-norm.txt";
		int[] pInds = new int[]{1, 2, 3};
		int[] fcInds = new int[]{10, 11, 12};

		convert(inFile, outFile, idIndex, pInds, fcInds);
	}

	private static void convertPDX() throws IOException
	{
		String inFile = "pdx.csv";
		int idIndex = 0;
		String outFile = "data-pdx.txt";
		int[] pInds = new int[]{1};
		int[] fcInds = new int[]{4};

		convert(inFile, outFile, idIndex, pInds, fcInds);
	}

	private static void convertPDXDrugT() throws IOException
	{
		String inFile = "pdx-drug-t.csv";
		int idIndex = 0;
		String outFile = "data-pdx-drug-t.txt";
		int[] pInds = new int[]{1, 2};
		int[] fcInds = new int[]{7, 8};

		convert(inFile, outFile, idIndex, pInds, fcInds);
	}

	private static void convertPDXDrug() throws IOException
	{
		String inFile = DATA_DIR + "PDX-drug.csv";
		int idIndex = 1;
		String outFile = BASE + "data-pdx-drug.txt";

		String[] h1 = FileUtil.readHeader(inFile, 2);
		String[] h2 = FileUtil.readHeader(inFile, 3);

		String[] header = new String[h1.length];
		for (int i = 0; i < header.length; i++)
		{
			header[i] = h1[i] + "_" + h2[i];
		}

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID");
		for (int i = 2; i < header.length; i++)
		{
			FileUtil.tab_write(header[i], writer);
		}

		FileUtil.linesTabbed(inFile).skip(4).forEach(t ->
		{
			FileUtil.lnwrite(t[idIndex].trim(), writer);
			for (int i = 2; i < t.length; i++)
			{
				FileUtil.tab_write(t[i], writer);
			}
		});

		writer.close();
	}

	private static void convertPerturbation() throws IOException
	{
		String inFile = "perturb.csv";
		int idIndex = 0;
		String outFile = "data-perturb.txt";
		int[] pInds = new int[]{1, 2, 3, 4};
		int[] fcInds = new int[]{13, 14, 15, 16};

		convert(inFile, outFile, idIndex, pInds, fcInds);
	}

	private static void convert(String inFile, String outFile, int idIndex, int[] pInds, int[] fcInds) throws IOException
	{
		inFile = DATA_DIR + inFile;
		String[] header = FileUtil.lines(inFile).findFirst().get().split("\t");

		outFile = BASE + outFile;
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID");
		for (int i = 0; i < pInds.length; i++)
		{
			writer.write("\t" + header[pInds[i]].replace("adj.P.Val.", ""));
		}

		FileUtil.linesTabbed(inFile).skip(1).forEach(t ->
		{
			FileUtil.lnwrite(t[idIndex].trim(), writer);
			for (int i = 0; i < pInds.length; i++)
			{
				String val = "";
				if (!t[pInds[i]].isEmpty())
				{
					val = Double.toString(Double.valueOf(t[pInds[i]]));

					if (t[fcInds[i]].startsWith("-")) val = "-" + val;
				}

				if (val.isEmpty()) val = "NaN";

				FileUtil.tab_write(val, writer);
			}
		});

		writer.close();

		Set<String> here = FileUtil.getTermsInTabDelimitedColumn(outFile, 0, 1);
		Set<String> inPlat = FileUtil.getTermsInTabDelimitedColumn(BASE + "platform.txt", 0, 1);
//		here = here.stream().map(s -> s.replaceAll(" ", "")).collect(Collectors.toSet());
//		inPlat = inPlat.stream().map(s -> s.replaceAll(" ", "")).collect(Collectors.toSet());
		CollectionUtil.printVennSets(here, inPlat);
	}

	static void prepareAnalysisFolders()
	{
		prepareAnalysisFolders("data-medullo.txt", "Medulloblastoma", "");
		prepareAnalysisFolders("data-pdx.txt", "PDX", "");
		prepareAnalysisFolders("data-perturb.txt", "Perturbation", "gene-activity = MAP2K1 i\ngene-activity = MAP2K2 i\ngene-activity = ALK i");
	}

	static void prepareAnalysisFolders(String dataFile, String rootDir, String extraParams)
	{
		String[] header = FileUtil.readHeader(BASE + dataFile);
		for (int i = 1; i < header.length; i++)
		{
			String dir = BASE + rootDir + "/" + header[i];
			FileUtil.createDirectories(dir);
			BufferedWriter writer = FileUtil.newBufferedWriter(dir + "/parameters.txt");
			FileUtil.write("proteomics-platform-file = ../../platform.txt\n" +
				"proteomics-values-file = ../../" + dataFile + "\n" +
				"id-column = ID\n" +
				"symbols-column = Symbols\n" +
				"sites-column = Sites\n" +
				"effect-column = Effect\n" +
				"\n" +
				"value-transformation = signed-p-values\n" +
				"\n" +
				"threshold-for-data-significance = 0.1 phosphoprotein\n" +
				"color-saturation-value = 10\n" +
				"show-all-genes-with-proteomic-data = true\n" +
				"\n" +
				"calculate-network-significance = true\n" +
				"use-network-significance-for-causal-reasoning = true\n" +
				"permutations-for-significance = 100000\n" +
				"fdr-threshold-for-network-significance = 0.1\n" +
				"\n" +
				"value-column = " + header[i] + "\n" + extraParams, writer);
			FileUtil.closeWriter(writer);
		}
	}

	private static void printMissingSites(String dir)
	{
		getMissingSiteEffectsRecursive(dir).stream().sorted().forEach(System.out::println);
	}

	private static Set<String> getMissingSiteEffectsRecursive(String parentDir)
	{
		if (Files.isDirectory(Paths.get(parentDir)))
		{
			Set<String> set = getMissingSiteEffects(parentDir);

			for (File dir : new File(parentDir).listFiles())
			{
				if (dir.isDirectory()) set.addAll(getMissingSiteEffectsRecursive(dir.getPath()));
			}

			return set;
		}
		return Collections.emptySet();
	}

	private static Set<String> getMissingSiteEffects(String dir)
	{
		System.out.println("dir = " + dir);
		Set<String> set = new HashSet<>();

		String sifFile = dir + "/causative.sif";
		String siteFile = dir + "/unknown-site-effects.txt";
		if (Files.exists(Paths.get(siteFile)) && Files.exists(Paths.get(sifFile)))
		{
			Set<String> regs = FileUtil.lines(sifFile).map(l -> l.split("\t")).filter(t -> t.length > 2).map(t -> t[0]).collect(Collectors.toSet());
			set.addAll(FileUtil.lines(siteFile).filter(s -> !regs.contains(s.substring(0, s.indexOf("_")))).collect(Collectors.toSet()));
		}

		return set;
	}
}
