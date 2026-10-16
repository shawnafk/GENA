import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def animate_gena_2d(fig, ax, coord_qq, coord_qq_minus_z, full_moon_pos, t, Ngap=60):
    # -------------------------- 1. Pre-slice data + validation --------------------------
    # Validate original data dimensions
    if coord_qq['GSM'].ndim != 2 or coord_qq['GSM'].shape[1] != 3:
        raise ValueError(f"coord_qq['GSM'] must be an N×3 array! Current shape: {coord_qq['GSM'].shape}")
    if coord_qq_minus_z['GSM'].ndim != 2 or coord_qq_minus_z['GSM'].shape[1] != 3:
        raise ValueError(f"coord_qq_minus_z['GSM'] must be an N×3 array! Current shape: {coord_qq_minus_z['GSM'].shape}")
    
    # Validate time series matches data length
    if len(t) != len(coord_qq['GSM']):
        raise ValueError(f"Time series length ({len(t)}) must match data length ({len(coord_qq['GSM'])})")
    
    # Pre-slice data and time series for global reuse
    coord_qq_slice = coord_qq['GSM'][::Ngap, :].copy()
    coord_minus_z_slice = coord_qq_minus_z['GSM'][::Ngap, :].copy()
    t_slice = t.iloc[::Ngap].copy()  # Slice time series to match animation frames
    
    # Validate sliced data
    num_frames = len(coord_qq_slice)
    if num_frames == 0:
        max_ngap = len(coord_qq['GSM']) - 1
        raise ValueError(f"Ngap={Ngap} is too large! Original data length is {len(coord_qq['GSM'])}, maximum valid Ngap is {max_ngap}")
    
    # -------------------------- 2. Element storage (closure variables) --------------------------
    static_artists = []  # Static elements (orbits, Earth, legend, etc.)
    dynamic_elements = []  # Dynamic elements (arrows and text)
    
    # -------------------------- 3. Initialize static elements --------------------------
    def init():
        # Clear previous static elements
        for artist in static_artists:
            artist.remove()
        static_artists.clear()
        
        # Plot GENA trajectory with pre-sliced data
        gena_line, = ax.plot(
            coord_qq_slice[:, 0], coord_qq_slice[:, 1],
            'r.', markersize=1.5, label='GENA'
        )
        
        # Plot Moon orbit and markers
        orbit_line, = ax.plot(
            full_moon_pos['GSM'][:, 0], full_moon_pos['GSM'][:, 1],
            'k', label='Moon Orbit'
        )
        marker1, = ax.plot(
            full_moon_pos['GSM'][::6, 0], full_moon_pos['GSM'][::6, 1],
            'k.', markersize=5
        )
        marker2, = ax.plot(
            full_moon_pos['GSM'][::24, 0], full_moon_pos['GSM'][::24, 1],
            'kx'
        )
        moon_init = ax.scatter(
            full_moon_pos['GSM'][0, 0], full_moon_pos['GSM'][0, 1],
            facecolors='k', s=100, label='Moon'
        )
        
        # Plot Earth
        earth = ax.scatter(
            0, 0, s=500, facecolors='none', edgecolors='b',
            label='Earth', linewidths=2
        )
        
        # Axis configuration
        ax.set_xlabel('GSM X (Re)')
        ax.set_ylabel('GSM Y (Re)')
        ax.set_title('GENA Satellite Orbit + Moon Orbit Animation')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.5)
        
        # Add legend
        legend = ax.legend(loc='upper right', framealpha=0.9)
        
        # Collect static elements
        static_artists.extend([gena_line, orbit_line, marker1, marker2, moon_init, earth, legend])
        return static_artists
    
    # -------------------------- 4. Update dynamic elements --------------------------
    def update(n):
        # Clear previous frame dynamic elements
        for elem in dynamic_elements:
            elem.remove()
        dynamic_elements.clear()
        
        # Get current frame data
        x, y = coord_qq_slice[n, 0], coord_qq_slice[n, 1]
        dx = coord_minus_z_slice[n, 0] * 10  # Scale direction vector for visibility
        dy = coord_minus_z_slice[n, 1] * 10
        current_time = t_slice.iloc[n]  # Current time for this frame
        
        # Draw arrows (add labels only for first frame)
        arrow1 = ax.arrow(
            x, y, -x, -y,
            head_width=0.2, head_length=0.1,
            fc='orange', ec='orange', linewidth=2,
            label='Towards Earth' if n == 0 else ""
        )
        arrow2 = ax.arrow(
            x, y, dx, dy,
            head_width=0.2, head_length=0.1,
            fc='green', ec='green', linewidth=2,
            label='Direction Vector' if n == 0 else ""
        )
        
        # Add text annotation for frame number and time
        time_text = ax.text(
            0.02, 0.98,  # Position (x=2%, y=98% in axes coordinates)
            f"frames = {n}, times = {str(current_time)}",  # Display text
            transform=ax.transAxes,  # Use axes coordinates (0-1 range)
            verticalalignment='top',  # Align text to top
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)  # White background with transparency
        )
        
        # Update legend for first frame to include arrow labels
        if n == 0:
            ax.legend(loc='upper right', framealpha=0.9)
        
        # Collect current frame dynamic elements
        dynamic_elements.extend([arrow1, arrow2, time_text])
        return dynamic_elements
    
    # -------------------------- 5. Create animation --------------------------
    ani = animation.FuncAnimation(
        fig=fig,
        func=update,
        frames=num_frames,
        init_func=init,
        blit=False,         # Set to False for better compatibility
        interval=100,       # Frame interval in milliseconds
        repeat=True         # Loop animation
    )
    
    return ani


# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 8))

# Call animation function with time series
ani = animate_gena_2d( fig, ax, coord_qq, coord_qq_minus_z, full_moon_pos, t=st_datetime, Ngap=20)

# Display animation
plt.show()
