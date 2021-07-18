package org.panda.misc.causalpath;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

public class CausalPathPenaltyBoxLister
{
	private static Set<String> getDisconnectedSymbols(String sifFile)
	{
		Set<String> covered = FileUtil.linesTabbed(sifFile).filter(t -> t.length > 2).map(t -> new String[]{t[0], t[2]})
			.flatMap(Arrays::stream).collect(Collectors.toSet());

		return FileUtil.linesTabbed(sifFile).filter(t -> t.length == 1).map(t -> t[0]).filter(s -> !covered.contains(s))
			.collect(Collectors.toSet());
	}

	private static Set<String> getIDs(Set<String> syms, String formatFile)
	{
		Set<String> ids = new HashSet<>();
		FileUtil.linesTabbed(formatFile).filter(t -> t[0].equals("node") && syms.contains(t[1])).forEach(t ->
		{
			if (t[2].equals("color")) ids.add(t[1]);
			else if (t[2].equals("rppasite"))
			{
				ids.add(t[3].substring(0, t[3].indexOf("|")));
			}
		});
		return ids;
	}

	private static Set<String> readValueLines(String valuesFile, Set<String> ids)
	{
		return FileUtil.lines(valuesFile).filter(l -> ids.contains(l.split("\t")[0])).collect(Collectors.toSet());
	}

	public static void writePenaltyBoxAsList(String dir) throws IOException
	{
		Set<String> disconnectedSymbols = getDisconnectedSymbols(dir + "/causative.sif");
		Set<String> ids = getIDs(disconnectedSymbols, dir + "/causative.format");
		Set<String> lines = readValueLines(dir + "/value-changes.txt", ids);
		BufferedWriter writer = FileUtil.newBufferedWriter(dir + "/penalty-box-values.txt");
		writer.write("ID\tChange\tAdjPval");
//		lines.stream().sorted().forEach(l -> FileUtil.lnwrite(l, writer));
		lines.stream().sorted().map(l -> l.split("\t")).forEach(t -> FileUtil.lnwrite(t[0] + "\t" + t[1] + "\t" + t[2], writer));
		writer.close();
	}

}
