package org.panda.misc.proteomics;

import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.statistics.Histogram;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by babur on 4/13/16.
 */
public class SIFRecurrenceCalculator
{
	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/RPPA/TCGA/correlation/";
//		prepareIntegratedSIFs(base, base + "recurrent");
		plotRecurrenceHistogram(base, base + "recurrent");
	}

	public static void prepareIntegratedSIFs(String inDir, String outDir) throws IOException
	{
		if (!(new File(outDir).exists())) new File(outDir).mkdirs();

		Map<String, Integer> sifCnt = new HashMap<>();
		Map<String, Integer> fmtCnt = new HashMap<>();
		Map<String, Integer> geneCnt = new HashMap<>();

		Arrays.stream(new File(inDir).listFiles())
			.filter(f -> f.isDirectory() && Files.exists(Paths.get(f.getPath() + "/causative.sif")))
			.map(file -> file.getPath() + "/causative")
			.forEach(pathWithoutExtension -> updateCounts(pathWithoutExtension, sifCnt, fmtCnt, geneCnt));


		int max = sifCnt.values().stream().reduce(0, Integer::max);
		int maxGene = geneCnt.values().stream().reduce(0, Integer::max);

		for (int i = max - 1; i > 0; i--)
		{
			ValToColor posVtc = new ValToColor(new double[]{i, max}, new Color[]{new Color(220, 255, 220), new Color(0, 80, 0)});
			ValToColor negVtc = new ValToColor(new double[]{i, max}, new Color[]{new Color(255, 220, 220), new Color(80, 0, 0)});
			ValToColor genVtc = new ValToColor(new double[]{i, maxGene}, new Color[]{new Color(220, 220, 220), new Color(0, 0, 0)});

			BufferedWriter w1 = new BufferedWriter(new FileWriter(outDir + "/" + i + ".sif"));

			final int thr = i;

			List<String> sifLines = sifCnt.keySet().stream().sorted(Comparator.comparing(sifCnt::get).reversed())
				.filter(line -> sifCnt.get(line) >= thr).collect(Collectors.toList());

			sifLines.forEach(line -> FileUtil.writeln(line, w1));

			w1.close();

			System.out.println("i = " + i);
			System.out.println("Phospho ratio = " + (sifLines.stream().filter(l -> l.contains("phospho")).count() / (double) sifCnt.keySet().stream().filter(l -> l.contains("phospho")).count()));
			System.out.println("Express ratio = " + (sifLines.stream().filter(l -> !l.contains("phospho")).count() / (double) sifCnt.keySet().stream().filter(l -> !l.contains("phospho")).count()));

			Files.createDirectories(Paths.get(outDir));

			BufferedWriter w2 = new BufferedWriter(new FileWriter(outDir + "/" + i + ".format"));
			w2.write("node\tall-nodes\tborderwidth\t2\n");
			w2.write("edge\tall-edges\twidth\t2\n");

			fmtCnt.keySet().stream().sorted(Comparator.comparing(fmtCnt::get).reversed())
				.filter(line -> fmtCnt.get(line) >= thr).forEach(line ->
				FileUtil.writeln(line, w2));

			Set<String> genes = sifLines.stream().map(l -> l.split("\t")).map(t -> new String[]{t[0], t[2]})
				.flatMap(Arrays::stream).collect(Collectors.toSet());

			genes.forEach(g -> FileUtil.writeln("node\t" + g + "\tbordercolor\t" + genVtc.getColorInString(geneCnt.get(g)), w2));
			sifLines.stream().forEach(l ->
			{
				String[] t = l.split("\t");

				FileUtil.writeln("edge\t" + t[0] + " " + t[1] + " " + t[2] + "\tcolor\t" +
					((t[1].startsWith("de") || t[1].startsWith("down") || t[1].startsWith("inh")) ?
						negVtc.getColorInString(sifCnt.get(l)) : posVtc.getColorInString(sifCnt.get(l))), w2);
			});

			w2.close();
		}
	}

	private static void updateCounts(String pathWithoutExtension, Map<String, Integer> sifCnt,
		Map<String, Integer> fmtCnt, Map<String, Integer> geneCnt)
	{
		Set<String> lines = FileUtil.lines(pathWithoutExtension + ".sif").collect(Collectors.toSet());
		lines.stream().forEach(line -> sifCnt.put(line, sifCnt.containsKey(line) ? sifCnt.get(line) + 1 : 1));

		Set<String> genes = lines.stream().map(l -> l.split("\t")).map(t -> new String[]{t[0], t[2]})
			.flatMap(Arrays::stream).collect(Collectors.toSet());

		genes.forEach(g -> geneCnt.put(g, geneCnt.containsKey(g) ? geneCnt.get(g) + 1 : 1));

		lines = FileUtil.lines(pathWithoutExtension + ".format").collect(Collectors.toSet());
		lines.stream().forEach(line -> fmtCnt.put(line, fmtCnt.containsKey(line) ? fmtCnt.get(line) + 1 : 1));
	}

	public static void plotRecurrenceHistogram(String inDir, String outDir) throws IOException
	{
		if (!(new File(outDir).exists())) new File(outDir).mkdirs();

		Map<String, Integer> sifCnt = new HashMap<>();
		Map<String, Integer> fmtCnt = new HashMap<>();
		Map<String, Integer> geneCnt = new HashMap<>();

		Arrays.stream(new File(inDir).listFiles())
			.filter(f -> f.isDirectory() && Files.exists(Paths.get(f.getPath() + "/causative.sif")))
			.map(file -> file.getPath() + "/causative")
			.forEach(pathWithoutExtension -> updateCounts(pathWithoutExtension, sifCnt, fmtCnt, geneCnt));

		int[] h = new int[sifCnt.values().stream().max(Integer::compareTo).get()];
		sifCnt.keySet().stream().filter(s -> !s.contains("phospho")).forEach(sif -> h[sifCnt.get(sif) - 1]++);
		for (int i = 0; i < h.length; i++)
		{
			System.out.println((i + 1) + "\t" + h[i]);
		}
	}
}