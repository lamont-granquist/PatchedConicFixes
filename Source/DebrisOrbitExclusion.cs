using HarmonyLib;
using UnityEngine;

// ReSharper disable ArrangeTypeModifiers
// ReSharper disable UnusedMember.Local

namespace PatchedConicFixes
{
    // Debris never needs to be considered for an encounter, a closest approach, or a
    // maneuver, but every consumer of Planetarium.Orbits still has to scan past it on
    // every call, since CheckEncounter and GetClosestApproach walk that whole list once
    // per orbit that shares a reference body with whatever is currently being solved.
    // In a save that has accumulated a lot of spent stages this adds up on its own, and
    // gets worse fast in a system with several moons where a lot of vessels and debris
    // end up sharing the same reference body.
    //
    // OrbitDriver.Start() is the only place anything is ever added to Planetarium.Orbits,
    // so pulling debris back out right after Start() finishes keeps every other consumer
    // of the list untouched, rather than trying to patch each call site separately.
    [HarmonyPatch(typeof(OrbitDriver), "Start")]
    class OrbitDriver__Start
    {
        static void Postfix(OrbitDriver __instance)
        {
            if (__instance.vessel != null && __instance.vessel.vesselType == VesselType.Debris)
            {
                Planetarium.Orbits.Remove(__instance);
            }
        }
    }
}
