# PatchedConicFixes

A Kerbal Space Program mod that fixes issues with patched conic trajectories.

## Requirements

- Kerbal Space Program 1.12.x

## Installation

Install via CKAN, or manually copy `GameData/PatchedConicFixes` into your KSP `GameData` folder.

## Building from Source

### Prerequisites

- [.NET SDK 6.0+](https://dotnet.microsoft.com/download)
- KSP 1.12.x installed

### Setup

Edit `Source/Directory.Build.props` to point to your KSP installation:

```xml
<KSPData>/path/to/KSP.app/Contents/Resources/Data</KSPData>     <!-- Mac -->
<KSPGameData>/path/to/KSP/GameData</KSPGameData>
```

Or create a `Source/Directory.Build.props.user` (gitignored) to keep your local paths out of version control:

```xml
<Project>
  <PropertyGroup>
    <KSPData>/your/local/path/KSP.app/Contents/Resources/Data</KSPData>
    <KSPGameData>/your/local/path/GameData</KSPGameData>
  </PropertyGroup>
</Project>
```

### Build

```
cd Source
dotnet build
```

The compiled DLL will be placed in `GameData/PatchedConicFixes/Plugins/`.

## License

MIT
