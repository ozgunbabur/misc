package org.panda.misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by babur on 4/13/16.
 */
public class SIFRecurrenceCalculator
{
	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/RPPA/TCGA/basic-correlation/";
		prepareIntegratedSIFs(base + "single", base + "integrated");
	}
	public static void prepareIntegratedSIFs(String inDir, String outDir) throws IOException
	{
		if (!(new File(outDir).exists())) new File(outDir).mkdirs();

		Map<String, Integer> sifCnt = new HashMap<>();
		Map<String, Integer> fmtCnt = new HashMap<>();

		Arrays.stream(new File(inDir).listFiles())
			.map(file -> file.getPath().substring(0, file.getPath().lastIndexOf(".")))
			.distinct().forEach(pathWithoutExtension -> updateCounts(pathWithoutExtension, sifCnt, fmtCnt));

		int max = sifCnt.values().stream().reduce(0, Integer::max);

		for (int i = max; i > 0; i--)
		{
			BufferedWriter w1 = new BufferedWriter(new FileWriter(outDir + "/" + i + ".sif"));

			final int thr = i;
			sifCnt.keySet().stream().sorted(Comparator.comparing(sifCnt::get).reversed())
				.filter(line -> sifCnt.get(line) >= thr).forEach(line -> {
				try
				{
					w1.write(line + "\n");
				} catch (IOException e)
				{
					throw new RuntimeException(e);
				}
			});

			w1.close();

			BufferedWriter w2 = new BufferedWriter(new FileWriter(outDir + "/" + i + ".format"));

			fmtCnt.keySet().stream().sorted(Comparator.comparing(fmtCnt::get).reversed())
				.filter(line -> fmtCnt.get(line) >= thr).forEach(line -> {
				try
				{
					w2.write(line + "\n");
				} catch (IOException e)
				{
					throw new RuntimeException(e);
				}
			});

			w2.close();
		}
	}

	private static void updateCounts(String pathWithoutExtension, Map<String, Integer> sifCnt,
		Map<String, Integer> fmtCnt)
	{
		try
		{
			Set<String> lines = Files.lines(Paths.get(pathWithoutExtension + ".sif")).collect(Collectors.toSet());
			lines.stream().forEach(line -> sifCnt.put(line, sifCnt.containsKey(line) ? sifCnt.get(line) + 1 : 1));
			lines = Files.lines(Paths.get(pathWithoutExtension + ".format")).collect(Collectors.toSet());
			lines.stream().forEach(line -> fmtCnt.put(line, fmtCnt.containsKey(line) ? fmtCnt.get(line) + 1 : 1));
		}
		catch (IOException e){ throw new RuntimeException(e);}
	}
}
