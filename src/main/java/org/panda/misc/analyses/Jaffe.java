package org.panda.misc.analyses;

import org.panda.causalpath.run.MergeResultsIntoSeriesView;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.statistics.UniformityChecker;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class Jaffe
{
	public static void main(String[] args) throws IOException
	{
//		generateDirectoriesRound2();
//		plotGeneVSEdgeCounts();
//		plotUniformityOfGraphSizePValues();
		evaluateNetworkSuccess();

//		CausalPathSubnetwork.printSignificantGenesRecursively("/home/ozgun/Analyses/Jaffe/round2/PC3", 0.1);

//		generateSeriesViews();
	}

	/**
	 * @deprecated keeping for historical reasons
	 */
	public static void generateDirectoriesRound1() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Jaffe/";
		String[] header = Files.lines(Paths.get(dir + "P29_CausalPath_Signed_31.txt")).findFirst().get().split("\t");

		for (int i = 7; i < header.length; i++)
		{
			String drug = header[i];

			FileUtil.mkdirs(dir + drug);
			BufferedWriter writer = FileUtil.newBufferedWriter(dir + drug + "/parameters.txt");

			Files.lines(Paths.get(dir + "parameters.txt.template")).forEach(l -> FileUtil.writeln(l, writer));

			FileUtil.writeln("value-column = " + drug, writer);

			writer.close();
		}
	}

	public static void generateDirectoriesRound2() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Jaffe/round2/";
		String cellType = "PC3";
//		String cellType = "A375";
//		String cellType = "MCF7";

		Map<String, List<String>> targets = Files.lines(Paths.get(dir + "drug-targets.txt")).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Arrays.asList(t[1].split(";"))));

		String[] header = Files.lines(Paths.get(dir + cellType + ".txt")).findFirst().get().split("\t");

		for (int i = 4; i < header.length; i++)
		{
			String drug = header[i];

			FileUtil.mkdirs(dir + cellType + "/" + drug);
			BufferedWriter writer = FileUtil.newBufferedWriter(dir + cellType + "/" + drug + "/parameters.txt");

			Files.lines(Paths.get(dir + "parameters.txt." + cellType + ".template")).forEach(l -> FileUtil.writeln(l, writer));

			for (String target : targets.getOrDefault(drug, Collections.emptyList()))
			{
				FileUtil.writeln("gene-activity = " + target + " i", writer);
			}

			FileUtil.lnwrite("value-column = " + drug, writer);

			writer.close();
		}
	}

	public static void plotGeneVSEdgeCounts() throws IOException
	{
		String base = "/home/ozgun/Analyses/Jaffe/round1";

		for (File dir : new File(base).listFiles())
		{
			if (dir.isDirectory())
			{
				int[] cnts = SIFFileUtil.getNodeAndEdgeCounts(dir + "/causative.sif");
				System.out.println(dir.getName() + "\t" + cnts[0] + "\t" + cnts[1]);
			}
		}
	}

	public static void plotUniformityOfGraphSizePValues()
	{
		String base = "/home/ozgun/Analyses/Jaffe/round2/PC3";
		ArrayList<Double> pvals = new ArrayList<>();
		collectP(new File(base), pvals);
		UniformityChecker.plot(pvals);
	}

	public static void collectP(File dir, List<Double> pvals)
	{
		String s = dir.getPath() + "/significance-pvals.txt";
		if (Files.exists(Paths.get(s)) &&
			(FileUtil.lines(dir.getPath() + "/causative.sif").findFirst().isPresent() &&
			FileUtil.lines(dir.getPath() + "/causative.sif").findFirst().get().contains("\t")))
		{
			String l = FileUtil.lines(s).findFirst().get();
			double p = Double.valueOf(l.substring(l.lastIndexOf(" ") + 1));
			pvals.add(p);
		}

		for (File sub : dir.listFiles())
		{
			if (sub.isDirectory()) collectP(sub, pvals);
		}
	}

	public static void evaluateNetworkSuccess() throws IOException
	{
		String root = "/home/ozgun/Analyses/CausalPath-paper/PC3";

		int[] cnts = {0, 0, 0};
		evaluateNetworkSuccessRecursive(new File(root), cnts);
		System.out.println("\nread = " + cnts[0]);
		System.out.println("pres = " + cnts[1]);
		System.out.println("enri = " + cnts[2]);
	}

	public static void evaluateNetworkSuccessRecursive(File dir, int[] cnts) throws IOException
	{
		if (Files.exists(Paths.get(dir.getPath() + "/parameters.txt")))
		{
			boolean[] b = checkIfGraphContainsInhibited(dir, cnts);
			System.out.println((b[0] ? "X" : " ") + "\t" + (b[1] ? "X" : " ") + "\t" + dir.getPath());
		}

		for (File sub : dir.listFiles())
		{
			if (sub.isDirectory()) evaluateNetworkSuccessRecursive(sub, cnts);
		}
	}



	private static boolean[] checkIfGraphContainsInhibited(File dir, int[] cnts) throws IOException
	{
		cnts[0]++;
		if (!Files.exists(Paths.get(dir.getPath() + "/results.txt"))) return new boolean[]{false, false};

		Set<String> inhGenes = Files.lines(Paths.get(dir.getPath() + "/parameters.txt"))
			.filter(l -> l.startsWith("gene-activity = ")).map(l -> l.split(" ")[2]).collect(Collectors.toSet());

		boolean b1 = Files.lines(Paths.get(dir.getPath() + "/results.txt")).map(l -> l.split("\t"))
			.anyMatch(t -> inhGenes.contains(t[0]) && t[4].endsWith("-inactive"));

		boolean b2 = Files.lines(Paths.get(dir.getPath() + "/causative.format")).map(l -> l.split("\t"))
			.anyMatch(t -> inhGenes.contains(t[1]) && t[2].equals("bordercolor") && t[3].equals("180 0 20"));

		if (b1) cnts[1]++;
		if (b2) cnts[2]++;

		return new boolean[]{b1, b2};
	}

	private static void generateSeriesViews() throws IOException
	{
		String dir = "/home/ozgun/Analyses/Jaffe/round2/";

		for (File subDir : new File(dir + "MCF7").listFiles())
		{
			String out = dir + "Combined/" + subDir.getName();
			Files.createDirectories(Paths.get(out));

			MergeResultsIntoSeriesView.run(out,
				subDir.getPath().replaceAll("MCF7", "A375"),
				subDir.getPath(),
				subDir.getPath().replaceAll("MCF7", "PC3"));
		}
	}
}