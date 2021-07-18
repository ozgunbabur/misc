package org.panda.misc.causalpath;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.util.Set;
import java.util.stream.Collectors;

public class CausalPathPlatformAndValuesFileComparison
{
	public static void main(String[] args)
	{
		printIDsIntersection();
	}

	private static void printIDsIntersection()
	{
		String dir = "/Users/ozgun/Documents/Analyses/SigPath-CausalPath/";
		String pFile = dir + "platform.txt";
		String vFile = dir + "data-medulo-mod-t.txt";
		String idCol = "ID";

		String[] pHeader = FileUtil.readHeader(pFile);
		String[] vHeader = FileUtil.readHeader(vFile);

		int pInd = ArrayUtil.indexOf(pHeader, idCol);
		int vInd = ArrayUtil.indexOf(vHeader, idCol);

		Set<String> pIDs = FileUtil.linesTabbedSkip1(pFile).map(t -> t[pInd]).collect(Collectors.toSet());
		Set<String> vIDs = FileUtil.linesTabbedSkip1(vFile).map(t -> t[vInd]).collect(Collectors.toSet());

		CollectionUtil.printNameMapping("Platform", "Values");
		CollectionUtil.printVennSets(pIDs, vIDs);
	}
}
