using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Emit;
using HarmonyLib;

namespace PatchedConicFixes
{
    // The stock interval uses true anomaly in [-pi, pi]. Subdivision first
    // unwraps a seam-crossing interval by adding 2*pi, so each interior point
    // must also wrap by a full turn. Subtracting pi moves it to the opposite
    // side of the orbit and can expand the child interval instead of shrinking it.
    [HarmonyPatch(typeof(Targeting.Interval), nameof(Targeting.Interval.Subdivide))]
    internal static class TargetingIntervalSubdivision
    {
        private static IEnumerable<CodeInstruction> Transpiler(IEnumerable<CodeInstruction> instructions)
        {
            List<CodeInstruction> code = instructions.ToList();
            var wrapConstants = new List<int>();
            for (int i = 0; i + 1 < code.Count; i++)
            {
                if (code[i].opcode == OpCodes.Ldc_R8 && code[i].operand is double value
                    && value == Math.PI && code[i + 1].opcode == OpCodes.Sub)
                    wrapConstants.Add(i);
            }

            // Validate the whole match before changing any instruction. A different
            // KSP implementation needs review rather than a partial numerical patch.
            if (wrapConstants.Count != 2)
                throw new InvalidOperationException("Targeting.Interval.Subdivide: expected two pi wrap subtractions; found "
                    + wrapConstants.Count + ". The installed KSP method requires review.");

            foreach (int index in wrapConstants)
                code[index].operand = 2 * Math.PI;
            return code;
        }
    }
}
