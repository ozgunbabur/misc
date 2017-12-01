package org.panda.misc;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class DifferenceOfSIFs
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/BigMech/Diff-view/";
		takeDifferenceUse3Columns(dir + "pc.sif", dir + "reach.sif");

//		String dir = "/home/babur/Projects/causalpath/";
//		printVennCounts3Columns(dir + "pc.sif", dir + "reach.sif");
	}

	public static void takeDifferenceBasic(String file1, String file2) throws IOException
	{
		String out1 = file1.substring(0, file1.lastIndexOf(".")) + "-diff-from-" + file2.substring(file2.lastIndexOf("/") + 1);
		String out2 = file2.substring(0, file2.lastIndexOf(".")) + "-diff-from-" + file1.substring(file1.lastIndexOf("/") + 1);
		String out3 = file1.substring(0, file1.lastIndexOf(".")) + "-and-" + file2.substring(file2.lastIndexOf("/") + 1);

		Set<String> set1 = Files.lines(Paths.get(file1)).collect(Collectors.toSet());
		Set<String> set2 = Files.lines(Paths.get(file2)).collect(Collectors.toSet());

		Set<String> setInt = new HashSet<>(set1);
		setInt.retainAll(set2);

		set1.removeAll(setInt);
		set2.removeAll(setInt);

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(out1));
		set1.forEach(l -> FileUtil.writeln(l, writer1));
		writer1.close();

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(out2));
		set2.forEach(l -> FileUtil.writeln(l, writer2));
		writer2.close();

		BufferedWriter writer3 = Files.newBufferedWriter(Paths.get(out3));
		setInt.forEach(l -> FileUtil.writeln(l, writer3));
		writer3.close();
	}

	public static void takeDifferenceUse3Columns(String file1, String file2) throws IOException
	{
		String out1 = file1.substring(0, file1.lastIndexOf(".")) + "-diff-from-" + file2.substring(file2.lastIndexOf("/") + 1);
		String out2 = file2.substring(0, file2.lastIndexOf(".")) + "-diff-from-" + file1.substring(file1.lastIndexOf("/") + 1);
		String out3 = file2.substring(0, file2.lastIndexOf(".")) + "-and-" + file1.substring(file1.lastIndexOf("/") + 1);

		Map<String, String> map1 = getRelToLineMap(file1);
		Map<String, String> map2 = getRelToLineMap(file2);

		Set<String> set1 = new HashSet<>(map1.keySet());
		Set<String> set2 = new HashSet<>(map2.keySet());
		Set<String> setInt = new HashSet<>(set1);
		setInt.retainAll(set2);

		set1.removeAll(setInt);
		set2.removeAll(setInt);

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(out1));
		set1.forEach(l -> FileUtil.writeln(map1.get(l), writer1));
		writer1.close();

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(out2));
		set2.forEach(l -> FileUtil.writeln(map2.get(l), writer2));
		writer2.close();

		BufferedWriter writer3 = Files.newBufferedWriter(Paths.get(out3));
		setInt.forEach(l -> FileUtil.writeln(map2.get(l), writer3));
		writer3.close();
	}

	private static Map<String, String> getRelToLineMap(String file1) throws IOException
	{
		return Files.lines(Paths.get(file1)).collect(Collectors.toMap(l ->
		{
			String[] t = l.split("\t");
			return t[0] + "\t" + t[1] + "\t" + t[2];
		}, l -> l));
	}

	public static void printVennCounts3Columns(String file1, String file2) throws IOException
	{
		CollectionUtil.printNameMapping(file1, file2);
		CollectionUtil.printVennCounts(
			Files.lines(Paths.get(file1)).collect(Collectors.toSet()),
			Files.lines(Paths.get(file2)).collect(Collectors.toSet())
		);
	}
}
