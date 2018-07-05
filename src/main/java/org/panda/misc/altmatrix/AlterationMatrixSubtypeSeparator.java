package org.panda.misc.altmatrix;

import org.panda.utility.FileUtil;

import javax.lang.model.type.TypeMirror;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class AlterationMatrixSubtypeSeparator
{
	public static void separate(String matrixFile, String outDir, Map<String, Set<String>> typeMap, int minAltThr) throws IOException
	{
		String[] header = Files.lines(Paths.get(matrixFile)).findFirst().get().split("\t");

		Map<String, BufferedWriter> writerMap = new HashMap<>();

		for (String type : typeMap.keySet())
		{
			String dir = outDir + File.separator + type;
			Files.createDirectories(Paths.get(dir));
			writerMap.put(type, Files.newBufferedWriter(Paths.get(dir + File.separator + "DataMatrix.txt")));
		}

		for (int i = 1; i < header.length; i++)
		{
			for (String type : typeMap.keySet())
			{
				if (typeMap.get(type).contains(header[i]))
				{
					writerMap.get(type).write("\t" + header[i]);
				}
			}
		}

		Files.lines(Paths.get(matrixFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String gene = t[0];
			Set<String> types = getTypesWhereTheGeneIsAltered(t, header, typeMap, minAltThr);

			types.forEach(type -> FileUtil.lnwrite(gene, writerMap.get(type)));

			for (int i = 1; i < t.length; i++)
			{
				for (String type : types)
				{
					if (typeMap.get(type).contains(header[i]))
					{
						FileUtil.tab_write(t[i], writerMap.get(type));
					}
				}
			}
		});

		writerMap.values().forEach(FileUtil::closeWriter);
	}

	private static Set<String> getTypesWhereTheGeneIsAltered(String[] data, String[] header,
		Map<String, Set<String>> typeMap, int minAltThr)
	{
		if (minAltThr < 1) return typeMap.keySet();

		Map<String, Integer> countMap = new HashMap<>();

		for (int i = 1; i < header.length; i++)
		{
			if (!data[i].equals("0"))
			{
				for (String type : typeMap.keySet())
				{
					if (typeMap.get(type).contains(header[i]))
					{
						countMap.put(type, countMap.containsKey(type) ? countMap.get(type) + 1 : 1);
					}
				}
			}
		}

		return countMap.keySet().stream().filter(type -> countMap.get(type) >= minAltThr).collect(Collectors.toSet());
	}

	public static void main(String[] args) throws IOException
	{
		// read types to samples
		Map<String, Set<String>> typeMap = Files.lines(Paths.get("/media/babur/6TB1/TCGA-pancan/sample-to-tissue-mapping.txt"))
			.map(l -> l.split("\t")).peek(t -> t[1] = t[1].contains("-") ? t[1] = t[1].substring(0, t[1].indexOf("-")) : t[1])
			.collect(Collectors.groupingBy(t -> t[1], Collectors.mapping(t -> t[0], Collectors.toSet())));

		String dir = "/media/babur/6TB1/TCGA-pancan/whole/";

		separate(dir + "DataMatrix.txt", dir + "types", typeMap, 2);
	}
}
