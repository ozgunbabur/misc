package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader_old;
import org.panda.utility.ValToColor;
import org.panda.utility.statistics.Summary;

import java.awt.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;

/**
 * Created by babur on 2/23/16.
 */
public class PanCanMutexIndividualSIFGenerator
{
	public static final double SCORE_THR = 0.01;
	public static final String OUTFILE = "results.sif";

	public static void main(String[] args) throws IOException
	{
		traverseRecursive(new File("/home/babur/Documents/PanCan/PanCan-results/"));
	}

	private static void traverseRecursive(File dir) throws IOException
	{
		if (new File(dir.getPath() + "/ranked-groups.txt").exists())
		{
			process(dir);
		}
		else
		{
			for (File file : dir.listFiles())
			{
				if (file.isDirectory()) traverseRecursive(file);
			}
		}
	}

	private static void process(File dir) throws IOException
	{
		Set<MutexReader.Group> groups = MutexReader.readMutexResults(dir.getPath());
		Map<String, Set<String>> map = new HashMap<>();
		groups.stream().filter(g -> g.score <= SCORE_THR).forEach(g ->
		{
			for (String g1 : g.genes)
			{
				for (String g2 : g.genes)
				{
					if (g1.compareTo(g2) < 0)
					{
						if (!map.containsKey(g1)) map.put(g1, new HashSet<>());
						map.get(g1).add(g2);
					}
				}
			}
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir.getPath() + "/" + OUTFILE));
		for (String source : map.keySet())
		{
			for (String target : map.get(source))
			{
				writer.write(source + "\tin-same-group\t" + target + "\n");
			}
		}
		writer.close();
	}
}
